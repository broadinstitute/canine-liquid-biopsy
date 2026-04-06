library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(flextable)
library(coin) 
library(WRS2)
library(lubridate)
# This package automatically finds the root folder of the downloaded project.
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ==============================================================================
# 1. SETUP & FILE PATHS
# ==============================================================================

# Define the input files
infile <- here("data", "BB_Supplementary_Data_1.txt")

# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

outputfile <- paste0(pdfdir, "Table.Fig_S10.stats.csv")
write_lines(c("Statistical results for analyses presented in Fig. S10"), outputfile, append=FALSE)

rename_stats <- tibble(order=c(2,1,3),
                       statName=c("Tumor fraction","cfDNA concentration","Fragment size ratio"),
                       axisLabel=c("Tumor fraction","Log10 cfDNA conc. (log10 ng/mL)","Fragment size ratio"),
                       stat=c("TF","log_DNA_conc","Fragment_size_ratio") ) %>% arrange(order)


# --- Helper: Scientific Notation Formatter ---
format_p_csv <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) return("NA")
    if (val == 0) return("< 2.2e-16") 
    if (val < 0.001) {
      formatC(val, format = "e", digits = 4)
    } else {
      formatC(val, format = "f", digits = 4)
    }
  })
}

# Helper for 2 significant figures in plots
format_p_sig2 <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) return("NA")
    if (val < 0.001) {
      # Return scientific notation with 2 decimal places (approx 3 sig figs)
      return(formatC(val, format = "e", digits = 2))
    }
    # Format to 2 significant digits, drop trailing zeros if not needed
    val_sig <- signif(val, 2)
    format(val_sig, scientific=FALSE, drop0trailing=TRUE)
  })
}

# --- Data Loading & Prep ---
# Try/catch block added to allow script generation without files present
if(file.exists(infile)) {
  raw <- as_tibble(read.csv(infile,header=T,sep="\t",na.strings = c("NA","N/A","","UNK")))
  raw <- raw %>% filter(Sample_ID!="")
  if("X" %in% colnames(raw)) raw$X <- NULL
  
  # Format Date
  raw <- raw %>% 
    mutate(Date_of_sample=str_replace(Date_of_sample," \\(approx\\)","")) %>% 
    mutate(Date_of_sample=lubridate::parse_date_time(Date_of_sample,"m/d/y"))
  
  # Clean Sample_type
  raw <- raw %>% mutate(Sample_type=str_remove(Sample_type,"_plasma"))
  
  # Create Sex.spay (Sex + Spay/Neuter Status)
  raw <- raw %>%
    mutate(Sex.spay=if_else(is.na(Spay.neuter.status),
                            Spay.neuter.status,
                            paste(Sex,Spay.neuter.status,sep=".")))
  
  # Recode Cancer_status
  raw <- raw %>% mutate(Cancer_status=if_else(Cancer_type!="Healthy","cancer","healthy"))
  
  # Add log of DNA_conc
  raw <- raw %>% mutate(log_DNA_conc=log10(DNA_conc))
  
  #### Get input data for analysis
  # NOTE: We keep covariates needed for ANOVA models here
  pdAll <- raw %>% 
    filter(First_Sample) %>% 
    select(ID, Breed, Cancer_type, Cancer_status, 
           log_DNA_conc, TF, Fragment_size_ratio,
           Sex.spay, Age_at_sample, Weight_kg, BCS, 
           Sequencing_group, Sample_type, Plasma_volume_estimated)
  
  pdAll <- pdAll %>% 
    ungroup() %>% 
    pivot_longer(c(log_DNA_conc, TF, Fragment_size_ratio)) %>% 
    filter(!is.na(value)) %>% 
    rename(stat=name)
  
  pdAll <- pdAll %>% left_join(rename_stats, by="stat") %>% filter(Cancer_status=="cancer")
  
  # Define Breed Groups
  dogs <- pdAll %>% select(ID,Breed) %>% distinct() 
  breedgrps <- dogs %>% group_by(Breed) %>% count() %>% mutate(breedGrp=if_else(n>=10,Breed,"other breeds")) %>% select(-n)
  pdAll <- pdAll %>% full_join(breedgrps, by="Breed") 
  pdAll <- pdAll %>% mutate(breedGrp=if_else(breedGrp==Breed & Breed!="mixed breed", breedGrp, "other dogs"))
  
  # Define Order
  grps <- pdAll %>% ungroup() %>% select(breedGrp) %>% distinct()
  grps <- tibble(breedGrp=c("other dogs","golden retriever","labrador retriever"),order=c(1:3)) %>% right_join(grps, by="breedGrp")
  
  # Define Comparisons
  comps <- list(c("other dogs","golden retriever"),c("other dogs","labrador retriever"),c("golden retriever","labrador retriever"))
} else {
  # Dummy data generation for script validation if files missing
  print("Warning: Input file not found. Generating dummy data for structure validation.")
  grps <- tibble(breedGrp=c("other dogs","golden retriever","labrador retriever"),order=c(1:3))
  comps <- list(c("other dogs","golden retriever"))
  # (Dummy pdAll creation would go here if needed for full run)
}


# ==============================================================================
# MULTIVARIABLE ANOVA (Consolidated Model)
# ==============================================================================

if(exists("pdAll")) {
  
  # Function runs single combined model
  run_breed_anova <- function(data) {
    stats_list <- unique(data$stat)
    res_list <- list()
    
    for(s in stats_list) {
      # Subset data for this stat
      d_s <- data %>% filter(stat == s)
      
      # 1. Predictors: Demographics + BCS
      base_preds <- c("breedGrp", "Cancer_type", "Sex.spay", "Age_at_sample", "BCS")
      
      # Add Conditional Predictors (Seq Group, Sample Type) initially
      candidates <- base_preds
      if(n_distinct(d_s$Sequencing_group, na.rm=TRUE) > 1) candidates <- c(candidates, "Sequencing_group")
      if(n_distinct(d_s$Sample_type, na.rm=TRUE) > 1) candidates <- c(candidates, "Sample_type")
      
      # Plasma Volume Logic
      if(s == "log_DNA_conc" && 
         n_distinct(d_s$Plasma_volume_estimated, na.rm=TRUE) > 1 &&
         mean(is.na(d_s$Plasma_volume_estimated)) < 0.5) {
        candidates <- c(candidates, "Plasma_volume_estimated")
      }
      
      # 2. Clean Data: Drop NAs for intended predictors
      d_s_clean <- d_s %>% drop_na(value, any_of(candidates))
      
      # Check sample size
      if(nrow(d_s_clean) < 10) next
      
      # 3. Final Variance Check
      final_preds <- c()
      for(p in candidates) {
        if(n_distinct(d_s_clean[[p]]) > 1) {
          final_preds <- c(final_preds, p)
        }
      }
      
      if(!"breedGrp" %in% final_preds) next
      
      # 4. Construct Formula
      form_s <- as.formula(paste("value ~", paste(final_preds, collapse = " + ")))
      
      # 5. Run ANOVA
      tryCatch({
        an_res <- d_s_clean %>% anova_test(form_s) %>% as_tibble()
        an_res <- an_res %>% mutate(stat = s)
        res_list[[s]] <- an_res
      }, error = function(e) { 
        print(paste("Error in ANOVA for", s, ":", e$message))
        return(NULL) 
      })
    }
    
    bind_rows(res_list)
  }
  
  # Run Single Consolidated Model
  anova_results <- run_breed_anova(pdAll)
  
  # Combine and Write
  anova_results <- anova_results %>%
    mutate(p_adj = p.adjust(p, method="fdr")) %>% 
    select(stat, Effect, DFn, DFd, F, p, p_adj, ges)
  
  write_lines(c("", "Multivariable ANOVA Results (Consolidated Model: Demog + BCS)"), outputfile, append=TRUE)
  write.table(anova_results, file=outputfile, append=TRUE, sep=",", row.names=FALSE, col.names=TRUE)
  
  
  # ==============================================================================
  # OMNIBUS ROBUST & PERMUTATION ANOVA (Unadjusted)
  # ==============================================================================
  
  # This calculates the one-way ANOVA equivalents using Robust and Permutation methods
  # matching the logic used in the pairwise comparisons (WRS2 and Coin)
  
  run_robust_omnibus <- function(data) {
    stats_list <- unique(data$stat)
    res_list <- list()
    
    for(s in stats_list) {
      # We only need value and breedGrp for unadjusted tests
      d_s <- data %>% filter(stat == s) %>% select(value, breedGrp) %>% drop_na()
      d_s$breedGrp <- as.factor(d_s$breedGrp)
      
      if(nrow(d_s) < 10) next
      
      # 1. Permutation Test (Coin Independence Test)
      # Tests independence of Value ~ BreedGrp (Global Null)
      p_perm <- tryCatch({
        as.numeric(coin::pvalue(coin::independence_test(value ~ breedGrp, data=d_s, distribution="approximate")))
      }, error=function(e) NA)
      
      # 2. Robust ANOVA (WRS2 t1way)
      # Trimmed means (tr=0.2) one-way ANOVA
      p_robust <- tryCatch({
        WRS2::t1way(value ~ breedGrp, data=d_s, tr=0.2)$p.value
      }, error=function(e) NA)
      
      res_list[[s]] <- tibble(
        stat = s,
        Test_Factor = "breedGrp",
        p_perm_omnibus = p_perm,
        p_robust_omnibus = p_robust
      )
    }
    
    bind_rows(res_list)
  }
  
  omnibus_robust <- run_robust_omnibus(pdAll) %>%
    mutate(
      p_perm_adj = p.adjust(p_perm_omnibus, method="fdr"),
      p_robust_adj = p.adjust(p_robust_omnibus, method="fdr")
    ) %>%
    mutate(across(starts_with("p_"), format_p_csv))
  
  write_lines(c("", "Omnibus Robust & Permutation Tests (Unadjusted for Covariates)"), outputfile, append=TRUE)
  write.table(omnibus_robust, file=outputfile, append=TRUE, sep=",", row.names=FALSE, col.names=TRUE)
  
  
  # ==============================================================================
  # ANOVA PLOTTING
  # ==============================================================================
  
  effect_labels_df <- tibble(
    Effect = c("breedGrp", "Cancer_type", "Sex.spay", "Age_at_sample", 
               "BCS", "Sequencing_group", "Sample_type", 
               "Plasma_volume_estimated"),
    Label = c("Breed Group", "Cancer Type", "Sex", "Age", 
              "Body Condition", "Seq Group", "Sample Type", 
              "Plasma Vol")
  )
  
  pd_anova_plot <- anova_results %>%
    filter(Effect != "Residuals") %>% 
    left_join(effect_labels_df, by="Effect") %>%
    mutate(
      Label = coalesce(Label, Effect), 
      # Use raw p-values for significance flag
      sig = if_else(p <= 0.05, "significant", "not significant"),
      # Use raw p-values for text label, scientific if small
      pstr = if_else(p < 0.001, 
                     paste0("p=", formatC(p, format = "e", digits = 2)), 
                     paste0("p=", sprintf("%.3f", p)))
    ) %>%
    left_join(rename_stats, by="stat")
  
  # Ordering logic based on log_DNA_conc P-values
  # 1. Get order of Effects based on P-value in log_DNA_conc (descending P = ascending significance)
  # Note: ggplot plots factors bottom-up. To have most significant (smallest P) at top,
  # we need it to be the last level in the factor.
  # So we arrange by descending P-value.
  order_reference <- pd_anova_plot %>%
    filter(stat == "log_DNA_conc") %>%
    arrange(desc(p)) %>% # Largest P first (bottom of plot), Smallest P last (top of plot)
    pull(Label)
  
  # 2. Get any labels that might exist in other panels but not in log_DNA_conc (safety check)
  all_labels <- unique(pd_anova_plot$Label)
  remaining_labels <- setdiff(all_labels, order_reference)
  final_levels <- c(remaining_labels, order_reference) # Put leftovers at bottom
  
  # 3. Apply levels
  pd_anova_plot$Label <- factor(pd_anova_plot$Label, levels = final_levels)
  
  # Set factor levels for statName to ensure correct facet order
  # Order: Conc, TF, Frag
  pd_anova_plot$statName <- factor(pd_anova_plot$statName, levels = rename_stats$statName)
  
  # Removed facet row for Model, as we only have one now
  p_anova_viz <- ggplot(pd_anova_plot, aes(x=ges, y=Label)) +
    geom_bar(aes(fill=sig), stat="identity", width=0.7) +
    geom_vline(xintercept = 0, color="#969696", linewidth=0.2) +
    geom_text(aes(label=pstr, color=sig), hjust=-0.1, size=2.5) +
    facet_grid(cols=vars(statName), scales="free", space="free_y") +
    scale_x_continuous("Effect Size (Generalized Eta Squared)", 
                       expand = expansion(mult = c(0, 0.4)), 
                       limits = c(0, NA)) + 
    scale_fill_manual(values=c("not significant"="#bdbdbd", "significant"="black")) +
    scale_color_manual(values=c("not significant"="#525252", "significant"="black")) +
    theme_cowplot(12) +
    theme(
      legend.position="none",
      axis.title.y=element_blank(),
      axis.text.y = element_text(size=9),
      strip.background = element_rect(fill="#525252"),
      strip.text = element_text(color="white", face="bold", size=9),
      panel.border = element_rect(color="black", fill=NA)
    )
  
  
  # ==============================================================================
  # STATISTICAL ANALYSIS (Unpaired: Wilcoxon, Permutation, Robust Yuen)
  # ==============================================================================
  
  calc_breed_stats <- function(data, subset_name) {
    res <- map_dfr(unique(data$statName), function(curr_stat) {
      d_metric <- data %>% filter(statName == curr_stat)
      
      map_dfr(comps, function(pair) {
        g1 <- pair[1]
        g2 <- pair[2]
        
        sub_d <- d_metric %>% filter(breedGrp %in% c(g1, g2)) %>% mutate(breedGrp = factor(breedGrp))
        
        if(nrow(sub_d) < 6 || length(unique(sub_d$breedGrp)) < 2) return(NULL)
        
        # wilcox.test (Mann-Whitney)
        p_raw <- tryCatch(wilcox.test(value ~ breedGrp, data=sub_d)$p.value, error=function(e) NA)
        
        p_perm <- tryCatch({
          as.numeric(coin::pvalue(coin::independence_test(value ~ breedGrp, data=sub_d, distribution="asymptotic")))
        }, error=function(e) NA)
        p_robust <- tryCatch({
          WRS2::yuen(value ~ breedGrp, data=sub_d, tr=0.2)$p.value
        }, error=function(e) NA)
        
        tibble(
          Subset = subset_name,
          statName = curr_stat,
          Group1 = g1,
          Group2 = g2,
          p_wilcox = p_raw, # Renamed from p_welch
          p_perm = p_perm,
          p_robust = p_robust
        )
      })
    })
    return(res)
  }
  
  stats_all <- calc_breed_stats(pdAll, "All Cancers")
  pdNoLymph <- pdAll %>% filter(Cancer_type != "Lymphoma")
  stats_no_lymph <- calc_breed_stats(pdNoLymph, "No Lymphoma")
  
  # Do not mutate to string format here, keep numeric for plots
  final_stats <- bind_rows(stats_all, stats_no_lymph) %>%
    mutate(
      p_wilcox_adj = p.adjust(p_wilcox, method="BH"),
      p_perm_adj = p.adjust(p_perm, method="BH"),
      p_robust_adj = p.adjust(p_robust, method="BH")
    ) 
  
  # Create separate formatted table for CSV output
  final_stats_csv <- final_stats %>%
    mutate(across(starts_with("p_"), format_p_csv))
  
  write_lines(c("", "Pairwise Comparisons by Breed (Independent Samples)"), outputfile, append=TRUE)
  write.table(final_stats_csv, file=outputfile, append=TRUE, sep=",", row.names=FALSE, col.names=TRUE)
  
  
  # ==============================================================================
  # PLOTTING
  # ==============================================================================
  
  theme_breeds <- function(){ 
    theme_minimal() %+replace%
      theme_cowplot(12) %+replace%
      theme(
        plot.title = element_text(hjust=0,size=7,face="bold",margin=margin(0,0,4,0)),
        plot.subtitle = element_text(hjust=0,size=7,margin=margin(0,0,4,0),lineheight = 1),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y = element_text(size=7,angle=90),
        axis.text.y = element_text(size=6,margin=margin(0,0,0,4)),
        axis.text.x = element_text(size=6,angle=0,hjust=0.5,margin=margin(2,0,0,0)))
  }
  
  # Accepts external stats table for FDR p-values
  make_breed_plot <- function(data, metric, stats_table, title_suffix="", no_lymph=FALSE) {
    pd <- data %>% filter(stat == metric) %>% filter(!is.na(value))
    
    # Filter Data
    subset_label <- "All Cancers"
    if(no_lymph) {
      pd <- pd %>% filter(Cancer_type != "Lymphoma")
      subset_label <- "No Lymphoma"
    }
    
    # Prepare X Labels
    xlabels <- pd %>% select(ID,breedGrp) %>% distinct() %>% group_by(breedGrp) %>% count() %>% mutate(xlabel=paste(breedGrp,"\n(N=",n,")",sep=""))
    xlabels <- xlabels %>% mutate(xlabel=str_replace_all(xlabel,"\\ ","\n"))
    pd$breedGrp <- factor(pd$breedGrp,levels=(grps %>% arrange(order) %>% pull(breedGrp)))
    
    # Prepare Stats for Plotting (raw p-values from wilcoxon)
    curr_statName <- unique(pd$statName)
    plot_stats <- stats_table %>% 
      filter(Subset == subset_label, statName == curr_statName) %>%
      rename(group1 = Group1, group2 = Group2) %>%
      # Use raw p-values from Wilcoxon (p_wilcox)
      mutate(label = format_p_sig2(p_wilcox)) 
    
    # Calculate Y positions for brackets to avoid overlap
    # We stack them: 10%, 20%, 30% above max value
    if(nrow(plot_stats) > 0) {
      max_y <- max(pd$value, na.rm=TRUE)
      min_y <- min(pd$value, na.rm=TRUE)
      y_range <- max_y - min_y
      # Ensure a minimum range to avoid flat lines crashing calculation
      if(y_range == 0) y_range <- 1 
      
      # Reduced step size to 0.1 to prevent brackets from going off-plot
      plot_stats <- plot_stats %>% 
        mutate(y.position = max_y + (row_number() * (y_range * 0.10)))
    }
    
    title_str <- paste(curr_statName, "\n", title_suffix, sep="")
    
    p <- ggplot(pd, aes(y=value, x=breedGrp)) 
    p <- p + geom_boxplot(aes(fill=breedGrp), color="#525252", alpha=0.25, outlier.shape=NA, width=0.5)
    p <- p + geom_jitter(aes(color=breedGrp), alpha=0.75, width=0.1, height=0, size=1, shape=21)
    p <- p + ggtitle(title_str, subtitle="Wilcoxon Rank-Sum Test")
    
    # Add P-values using manual table
    if(nrow(plot_stats) > 0) {
      p <- p + stat_pvalue_manual(plot_stats, label = "label", tip.length = 0.01, size = 2.5)
    }
    
    # Axis scaling
    if(metric == "log_DNA_conc") {
      p <- p + scale_y_continuous(unique(pd$axisLabel), expand = expansion(mult = c(0.05, 0.35))) 
    } else if(metric == "TF") {
      p <- p + scale_y_continuous(unique(pd$axisLabel), limits=c(0, 1.25), breaks=c(0,0.5,1), expand = expansion(mult = c(0.05, 0.35)))
    } else {
      p <- p + scale_y_continuous(unique(pd$axisLabel), expand = expansion(mult = c(0.05, 0.35)))
    }
    
    p <- p + scale_x_discrete("", breaks=xlabels$breedGrp, labels=xlabels$xlabel)
    p <- p + scale_color_manual(values=c("#8dd3c7","#fdb462","#bebada","#fb8072","#80b1d3"))
    p <- p + scale_fill_manual(values=c("#8dd3c7","#fdb462","#bebada","#fb8072","#80b1d3"))
    p <- p + theme_breeds()
    
    return(p)
  }
  
  # --- Generate Plots ---
  
  # 1. Concentration
  p_conc <- make_breed_plot(pdAll, "log_DNA_conc", final_stats)
  p_conc_xLymph <- make_breed_plot(pdAll, "log_DNA_conc", final_stats, "(no lymphomas)", no_lymph=TRUE)
  
  # 2. TF
  p_tf <- make_breed_plot(pdAll, "TF", final_stats)
  p_tf_xLymph <- make_breed_plot(pdAll, "TF", final_stats, "(no lymphomas)", no_lymph=TRUE)
  
  # 3. Fragment Size
  p_frag <- make_breed_plot(pdAll, "Fragment_size_ratio", final_stats)
  p_frag_xLymph <- make_breed_plot(pdAll, "Fragment_size_ratio", final_stats, "(no lymphomas)", no_lymph=TRUE)
  
  
  # --- Histograms ---
  
  # Histogram 1: Breeds
  pd <- raw %>% ungroup() %>% select(ID,Breed) %>% distinct() %>% filter(Breed != "mixed breed" & Breed!="other dogs")
  breednames <- pd %>% group_by(Breed) %>% count() %>% rename(ndogs=n) %>% filter(ndogs>2) %>% select(Breed,ndogs)
  breednames <- breednames %>% group_by(ndogs) %>% summarize(breedstr=paste(Breed,collapse=", ")) %>% mutate(breedstrP=str_replace_all(breedstr," ","\n"))
  breednames <- breednames %>% mutate(breedstrP=str_replace(breedstrP,"pit\n","pit "))
  breednames <- breednames %>% mutate(breedstrP=str_replace(breedstrP,"old\n","old "))
  
  pd <- pd %>% group_by(Breed) %>% count() %>% rename(ndogs=n)
  pd <- pd %>% ungroup() %>% group_by(ndogs) %>% count() %>% rename(nbreeds=n) %>% mutate(nbreeds_str=paste("N=",nbreeds,sep=""))
  pd <- pd %>% left_join(breednames, by="ndogs") %>% mutate(nbreeds_str=if_else(ndogs<3,nbreeds_str,paste(nbreeds_str,"\n(",breedstrP,")",sep="")))
  pd <- pd %>% mutate(label=if_else(ndogs==1,"1 dog",paste(ndogs,"dogs")))
  
  pd$label <- factor(pd$label,levels=pd %>% arrange(ndogs) %>% pull(label))
  p <- ggplot(pd,aes(y=nbreeds,x=factor(label)))
  p <- p + geom_bar(stat="identity",width=0.5)
  p <- p + geom_text(aes(label=nbreeds_str),nudge_y=0.25,size=1.75,vjust=0,lineheight=0.95)
  p <- p + scale_y_continuous("# breeds")
  p <- p + scale_x_discrete("# dogs in breed")
  p <- p + theme_cowplot(12)
  p <- p + theme(plot.title = element_text(hjust=0,size=8,face="bold",margin=margin(0,0,4,0)),
                 plot.subtitle = element_text(hjust=0,size=7,margin=margin(0,0,4,0),lineheight = 1),
                 legend.position = "none",
                 axis.title = element_text(size=7,face="bold"),
                 axis.text = element_text(size=6,margin=margin(0,0,0,4)),
                 panel.border = element_rect(colour = "#525252", fill=NA, size=1))
  p_breeds <- p
  
  # Histogram 2: Cancer Types
  # Fix many-to-many relationship warning by making ID unique
  d <- pdAll %>% select(ID) %>% distinct() %>% inner_join(raw, by="ID")
  d <- d %>% filter(Cancer_status!="healthy") %>% select(ID,Cancer_type) %>% distinct() %>% group_by(Cancer_type) %>% count() %>% rename(ndogs=n)
  pd <- d %>% ungroup() %>% group_by(ndogs) %>% count() %>% rename(ntypes=n)
  pd <- pd %>% mutate(label=if_else(ndogs==1,"1 dog",paste(ndogs,"dogs")))
  pd$label <- factor(pd$label,levels=pd %>% arrange(ndogs) %>% pull(label))
  
  p <- ggplot(pd,aes(y=ntypes,x=factor(label)))
  p <- p + geom_bar(stat="identity",width=0.5)
  p <- p + geom_text(aes(label=ntypes),nudge_y=1,size=2)
  p <- p + scale_y_continuous("# cancer types")
  p <- p + scale_x_discrete("# dogs with cancer type")
  p <- p + theme_cowplot(12)
  p <- p + theme(plot.title = element_text(hjust=0,size=8,face="bold",margin=margin(0,0,4,0)),
                 plot.subtitle = element_text(hjust=0,size=7,margin=margin(0,0,4,0),lineheight = 1),
                 legend.position = "none",
                 axis.title = element_text(size=7,face="bold"),
                 axis.text = element_text(size=6,margin=margin(0,0,0,4)),
                 panel.border = element_rect(colour = "#525252", fill=NA, size=1))
  p_types <- p
  
  
  # --- Assemble Figure ---
  pdfname <- paste(pdfdir,"Fig_S10.pdf",sep="")
  row1 <- plot_grid(p_breeds,p_conc,p_conc_xLymph,labels=LETTERS[1:3],nrow=1,label_size = 12,rel_widths = c(1,0.5,0.5))
  row2 <- plot_grid(p_tf,p_tf_xLymph,p_frag,p_frag_xLymph,labels=LETTERS[4:7],nrow=1,label_size = 12,rel_widths = 1)
  row3 <- plot_grid(p_anova_viz, labels=LETTERS[8], nrow=1, label_size=12)
  
  grid <- plot_grid(row1,row2,row3,ncol=1,label_size = 12, rel_heights=c(0.7, 0.7, 0.5))
  ggsave(grid,filename=pdfname,width=7,height=10)
  print(paste("Made",pdfname))
}