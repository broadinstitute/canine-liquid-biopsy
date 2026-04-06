library(ggpubr)
library(tidyverse)
library(rstatix)
library(cowplot)
library(flextable)
library(coin) 
library(WRS2) 
# This package automatically finds the root folder of the downloaded project.
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ==============================================================================
# 1. SETUP & FILE PATHS
# ==============================================================================

# Define the input files
infile <- here("data", "BB_Supplementary_Data_1.txt")
infile_human <- here("data", "HUMAN_DATA.txt")

# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

outputfile <- paste0(pdfdir, "Table.Fig_S3.stats.csv")
write_lines(c("Statistical results for analyses presented in Fig. S3"), outputfile, append=FALSE)

# Setup
rename_stats <- tibble(
  order = c(2, 1, 3),
  statName = c("Tumor fraction", "cfDNA concentration", "Fragment size ratio"),
  axisLabel = c("Tumor fraction", "cfDNA concentration (ng/mL)", "Fragment size ratio"),
  stat = c("TF", "log_DNA_conc", "Fragment_size_ratio")
) %>% arrange(order)

# --- Helper: Scientific Notation Formatter for CSV ---
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

# --- Helper: Permutation p for unpaired difference in means (Permuted T-Test) ---
.perm_p_indep <- function(x, g, B = 10000, seed = 1) {
  set.seed(seed)
  g <- as.factor(g)
  if (nlevels(g) != 2) return(NA)
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(g[ok])
  if(length(unique(g)) < 2) return(NA)
  
  d_obs <- diff(tapply(x, g, mean, na.rm = TRUE)) # group2 - group1
  
  perm_stats <- replicate(B, {
    gp <- sample(g)
    # Check if split is valid
    means <- tapply(x, gp, mean, na.rm = TRUE)
    if(length(means) < 2) return(FALSE)
    d  <- diff(means)
    abs(d) >= abs(d_obs)
  })
  (sum(perm_stats) + 1) / (B + 1)
}

# --- Data Loading & Cleaning ---
# get human data
human_pdAll <- as_tibble(read.csv(infile_human, header=T, na.strings=c("NA","N/A",""))) %>% filter(!is.na(First_Sample))
human_pdAll <- human_pdAll %>% mutate(value=if_else(stat=="DNA_conc",log10(value),value)) %>% mutate(stat=if_else(stat=="DNA_conc","log_DNA_conc",stat)) 
human_pdAll <- human_pdAll %>% mutate(Cancer_status2=Cancer_status, species="human")
human_pdAll <- human_pdAll %>% mutate(Sample_type="EDTA_plasma")

# get dog data
# Parsing safety - added quote="" and comment.char="" to prevent shifts, strip.white to handle spaces
dog_pdAll <- as_tibble(read.csv(infile, header=T, sep="\t", na.strings = c("NA","N/A",""), quote = "", comment.char = "", strip.white = TRUE)) 

# Clean Cancer_type
dog_pdAll <- dog_pdAll %>% 
  mutate(Cancer_type = trimws(Cancer_type))

# Logic to define Cancer_status 
dog_pdAll <- dog_pdAll %>% mutate(Cancer_status=if_else(Cancer_type!="Healthy","cancer","healthy"))

dog_pdAll <- dog_pdAll %>% mutate(Cancer_status2=Cancer_status)
dog_pdAll <- dog_pdAll %>% mutate(log_DNA_conc=log10(DNA_conc)) %>% select(-DNA_conc)
dog_pdAll <- dog_pdAll %>% mutate(species="dog")

# Allow "High" replicates as valid data points
# "NA" string check
dog_pdAll <- dog_pdAll %>% 
  filter(is.na(replicate) | replicate %in% c("Rep1", "High", "NA")) %>% 
  select(ID,species,Cancer_status,Cancer_type,Plasma_volume_estimated,First_Sample,Sample_type,log_DNA_conc,TF,Fragment_size_ratio)

dog_pdAll <- dog_pdAll %>% ungroup()  %>% pivot_longer(c(log_DNA_conc,TF,Fragment_size_ratio)) %>% filter(!is.na(value)) %>% rename(stat=name)

# Combine human and dog data
# Keep filter(First_Sample)
pdAll <- dog_pdAll %>% 
  bind_rows(human_pdAll %>% mutate(Plasma_volume_estimated=FALSE) %>% select(colnames(dog_pdAll))) %>% 
  filter(value!=-Inf) %>%
  left_join(rename_stats) %>%
  mutate(xpos=if_else(Cancer_status=="healthy",3,if_else(Cancer_status=="cancer",1,0))) %>%
  mutate(xpos=if_else(species=="human",xpos-0.5,xpos)) %>%
  # --- GLOBAL FILTER: FIRST SAMPLES ONLY ---
  filter(First_Sample)


# ==============================================================================
# STATISTICAL ANALYSIS BLOCK
# ==============================================================================

# --- Panel A: ANOVA + permuted p + robust p (dogs + humans) ---
letter_i <- 1
write_lines(
  c("", paste(LETTERS[letter_i], "ANOVA analysis on liquid biopsy metrics by species & cancer type (first samples only)")),
  outputfile,
  append = TRUE
)
letter_i <- letter_i + 1

# 1. Prepare Data: Collapse rare cancer types
min_N <- 10

pdAnova <- pdAll %>%
  filter(Cancer_type != "Unknown") %>%
  group_by(species, Cancer_type) %>%
  mutate(n_group = n()) %>%
  ungroup() %>%
  mutate(Cancer_type = if_else(n_group < min_N, "Other", Cancer_type)) %>%
  select(-n_group)

# 2. Define valid_groups for later use (for plots)
valid_groups <- pdAnova %>% 
  select(species, Cancer_type) %>% 
  distinct()

# 3. Run Statistical Tests
anova_results <- pdAnova %>%
  group_by(statName) %>%
  group_modify(~ {
    statname_curr <- unique(.y$statName)
    
    # --- Standard ANOVA ---
    base_res <- tryCatch({
      rstatix::anova_test(.x, value ~ Cancer_type + species) %>% as_tibble()
    }, error = function(e) {
      return(tibble(
        Effect = NA, DFn = NA, DFd = NA, F = NA, p = NA, `p<.05` = NA, ges = NA
      ))
    })
    
    if (nrow(base_res) == 0 || all(is.na(base_res$p))) {
      return(base_res %>% mutate(p_perm = NA, p_robust = NA))
    }
    
    # Setup for Permutation & Robust Tests
    .x <- .x %>%
      mutate(
        Cancer_type = factor(Cancer_type),
        species = factor(species),
        interaction_group = interaction(species, Cancer_type, drop = TRUE)
      )
    
    # --- Permuted p-value ---
    perm_p <- tryCatch({
      perm_mod <- coin::independence_test(
        value ~ interaction_group,
        data = .x,
        distribution = "asymptotic" 
      )
      as.numeric(coin::pvalue(perm_mod))
    }, error = function(e) {
      NA_real_
    })
    
    # --- Robust ANOVA (Trimmed Means) ---
    robust_p <- tryCatch({
      if (length(unique(.x$interaction_group)) < 2) {
        NA_real_
      } else {
        WRS2::t1way(value ~ interaction_group, data = .x, tr = 0.2)$p.value
      }
    }, error = function(e) {
      NA_real_
    })
    
    # Combine results
    base_res %>%
      mutate(
        p_perm = perm_p,
        flag_perm_disagree = ifelse(!is.na(p) & !is.na(perm_p), (p < 0.05) != (perm_p < 0.05), NA),
        p_robust = robust_p,
        flag_robust_disagree = ifelse(!is.na(p) & !is.na(robust_p), (p < 0.05) != (robust_p < 0.05), NA)
      )
  }) %>%
  ungroup()

# Apply formatting to CSV output only
anova_results_csv <- anova_results %>%
  mutate(
    p = format_p_csv(p),
    p_perm = format_p_csv(p_perm),
    p_robust = format_p_csv(p_robust)
  )

write.table(anova_results_csv, file = outputfile, append = TRUE, sep = ",", row.names = FALSE)


# ==============================================================================
# PANEL B-D STATS (Specific Types Only)
# ==============================================================================

write_lines(
  c("", paste(LETTERS[letter_i], "Pairwise comparisons (Specific Types Only) - Wilcoxon & Permuted T-test")),
  outputfile,
  append = TRUE
)
letter_i <- letter_i + 1

# Iterate over unique stats
stats_list <- unique(pdAll$statName)

pairwise_results <- map_dfr(stats_list, function(curr_stat) {
  
  # 1. Filter data for this stat
  d_stat <- pdAll %>%
    filter(statName == curr_stat, !is.na(value))
  
  # 2. Identify Valid Groups (N >= min_N OR Healthy)
  valid_types <- d_stat %>%
    filter(Cancer_type != "Unknown") %>%
    group_by(species, Cancer_type) %>%
    summarize(n = n(), .groups = "drop") %>%
    filter(n >= min_N | Cancer_type == "Healthy")
  
  # 3. Filter data to these valid types only
  d_stat <- d_stat %>%
    semi_join(valid_types, by = c("species", "Cancer_type")) %>%
    mutate(group = paste(species, Cancer_type, sep = " : "))
  
  if(length(unique(d_stat$group)) < 2) return(NULL)
  
  # 4. Generate combinations
  combos <- combn(unique(d_stat$group), 2, simplify = FALSE)
  
  # 5. Run tests
  map_dfr(combos, function(g) {
    g1 <- g[1]
    g2 <- g[2]
    sub_d <- d_stat %>% filter(group %in% c(g1, g2))
    sub_d$group <- factor(sub_d$group)
    
    # Wilcoxon Rank Sum Test
    p_wilcox <- tryCatch(wilcox.test(value ~ group, data=sub_d)$p.value, error=function(e) NA)
    
    # Permutation T-test (Difference in means)
    p_perm_t <- tryCatch({
      .perm_p_indep(sub_d$value, sub_d$group, B = 10000)
    }, error=function(e) NA)
    
    tibble(
      statName = curr_stat,
      group1 = g1,
      group2 = g2,
      wilcox_p = p_wilcox,
      perm_t_p = p_perm_t
    )
  })
})

if(nrow(pairwise_results) > 0){
  # FDR correction (grouped by statName)
  pairwise_results <- pairwise_results %>%
    group_by(statName) %>%
    mutate(wilcox_fdr = p.adjust(wilcox_p, method = "fdr")) %>%
    ungroup()
  
  pairwise_results_csv <- pairwise_results %>%
    mutate(
      wilcox_p = format_p_csv(wilcox_p),
      wilcox_fdr = format_p_csv(wilcox_fdr),
      perm_t_p = format_p_csv(perm_t_p)
    )
  write.table(pairwise_results_csv, file = outputfile, append = TRUE, sep = ",", row.names = FALSE)
} else {
  write_lines("No valid pairwise comparisons found.", outputfile, append = TRUE)
}


# ==============================================================================
# PLOTTING (Labels still include "Other")
# ==============================================================================

# --- Make Panel A Table (Cleaned for PDF) ---
format_p_display <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) return("NA")
    if (val < 0.001) {
      formatC(val, format = "e", digits = 2)
    } else {
      formatC(val, format = "f", digits = 3)
    }
  })
}

ft_data <- anova_results %>%
  select(statName, Effect, DFn, DFd, F, p, ges) %>%
  mutate(p = format_p_display(p))

ft <- flextable(ft_data)
ft <- ft %>% 
  theme_vanilla() %>% 
  autofit() %>% 
  colformat_double(j=c("F", "ges"), digits = 2) 

ft_plot <- gen_grob(ft, fit = "width", scaling="full", just = "center", width=0.2)


plot_boxes <- function(titlestr, statIn, cancerlabels, no_estimated){
  
  pd <- pdAll %>% filter(stat==statIn) %>% filter(!is.na(value)) 
  
  if (no_estimated){
    pd <- pd %>% filter(!Plasma_volume_estimated)
  }
  
  plotdata <- pd
  plotdata <- plotdata %>% left_join(cancerlabels) # Join mapping
  
  # Handle rare cancers that aren't in 'cancerlabels' mapping
  # Map them to "Other (Rare)".
  plotdata <- plotdata %>% 
    mutate(Cancer_type = if_else(is.na(new), "Other (Rare)", new)) %>% 
    select(-new) %>% 
    mutate(Cancer_status2 = paste(species, Cancer_type))
  
  xpositions <- plotdata %>% group_by(species,Cancer_status2) %>% summarize(ntotal=n(),mean=mean(value)) %>% arrange(mean) %>% ungroup() %>% mutate(xpos=row_number()) %>% select(-mean)
  plotdata <- plotdata %>% select(-xpos) %>% right_join(xpositions)
  xlabels <- plotdata %>% group_by(species,Cancer_status,Cancer_status2,xpos) %>% count() %>% mutate(xlabel=paste(tolower(str_replace(Cancer_status2," \\(","\n(")),"\nN=",n,sep="")) %>% select(-n)
  xlabels <- xlabels %>% mutate(xlabel=str_replace_all(xlabel,"\\ ","\n"))
  xlabels <- xlabels %>% mutate(xlabel=str_replace_all(xlabel,"\\_"," "))
  xlabels <- xlabels %>% mutate(xlabel=str_replace_all(xlabel,"osarcoma","o-\nsarcoma"))
  xlabels <- xlabels %>% mutate(xlabel=if_else(str_detect(xlabel,"metastatic"),xlabel,paste(xlabel,"\n")))
  
  Cancer_status2_levels <- plotdata %>% select(Cancer_status2,xpos) %>% distinct() %>% arrange(xpos) %>% pull(Cancer_status2)
  plotdata$Cancer_status2 <- factor(plotdata$Cancer_status2,levels=Cancer_status2_levels)
  
  statstrings <- plotdata %>% group_by(Cancer_status,Cancer_status2,species) %>% summarize(mean=mean(value),sd=sd(value)) %>% mutate(statstr=paste(round(mean,3),"\n+/-",round(sd,3),sep=""))
  
  rangey <- max(plotdata$value)-min(plotdata$value)
  miny <- min(plotdata$value) -  (rangey/20)
  
  p <- ggplot(plotdata,aes(y=value,x=Cancer_status2)) 
  p <- p + geom_text(aes(label=statstr,color=species),y=miny, data=statstrings,size=1.75,lineheight=0.9,hjust=0.5,vjust=1)
  p <- p + geom_boxplot(aes(fill=species),color="#525252",alpha=0.25,outlier.shape=16,outlier.alpha=0.5,outlier.size=1,width=0.5,linewidth=0.5)
  
  p <- p + ggtitle(unique(plotdata$statName),subtitle=titlestr)
  if (statIn=="log_DNA_conc"){
    p <- p + scale_y_continuous(unique(plotdata$axisLabel),breaks=c(-1:3),labels=10**c(-1:3),expand = expansion(mult = c(0.15, 0.1)))
  } else {
    p <- p + scale_y_continuous(unique(plotdata$axisLabel),expand = expansion(mult = c(0.15, 0.1)))
  }
  p <- p + scale_x_discrete("",breaks=xlabels$Cancer_status2,labels=xlabels$xlabel) 
  p <- p + scale_shape_manual(values=c(21,19))
  p <- p + scale_color_manual(values=c("#cb181d","#7570b3"))
  p <- p + scale_fill_manual(values=c("#cb181d","#7570b3"))
  p <- p + theme_minimal() %+replace%
    theme_cowplot(12) %+replace%
    theme(
      plot.title = element_text(hjust=0,size=7,face="bold",margin=margin(0,0,4,0)),
      plot.subtitle = element_text(hjust=0,size=6,margin=margin(0,0,4,0),lineheight = 1),
      legend.position = "none",
      axis.title.x=element_blank(),
      axis.title.y = element_text(size=6,angle=90),
      axis.text.y = element_text(size=5,margin=margin(0,0,0,4)),
      axis.text.x = element_text(size=5.5,angle=0,hjust=0.5,margin=margin(2,0,0,0),vjust=1))
}

# Generate Boxplot Mapping Using ALL DATA
# This ensures even rare cancers (which are excluded from stats via valid_groups) 
# get mapped to "Other"
type_counts <- pdAll %>%
  filter(Cancer_type != "Unknown") %>%
  group_by(stat, species, Cancer_type) %>%
  count(name = "ntot") %>%
  ungroup()

type_mapping <- type_counts %>%
  mutate(new_label_base = if_else(ntot >= min_N | Cancer_type == "Healthy", Cancer_type, "Other"))

group_counts <- type_mapping %>%
  group_by(stat, species, new_label_base) %>%
  summarize(ntypes = n(), .groups = "drop")

cancertypes <- type_mapping %>%
  left_join(group_counts, by = c("stat", "species", "new_label_base")) %>%
  mutate(new = if_else(new_label_base == "Other", 
                       paste0("Other (", ntypes, "_types)"), 
                       new_label_base)) %>%
  select(stat, species, Cancer_type, new)

# Plot panels
# no_estimated = FALSE (Keep samples with estimated volumes)
conc_boxes <- plot_boxes("First samples", "log_DNA_conc", cancertypes, no_estimated = FALSE)
tf_boxes   <- plot_boxes("First samples", "TF", cancertypes, no_estimated = FALSE)
frag_boxes <- plot_boxes("First samples", "Fragment_size_ratio", cancertypes, no_estimated = FALSE)


# --- Volume vs Concentration ---
pd <- pdAll %>% filter(stat=="log_DNA_conc") %>% filter(!is.na(value)) %>% filter(species=="dog")
pd <- pd %>% mutate(Cancer_type="All") %>% bind_rows(pd)
bad_types <- pd %>% group_by(Cancer_type,Plasma_volume_estimated) %>% count() %>% ungroup() %>% filter(n<=5) %>% select(Cancer_type) %>% distinct()
pd <- pd %>% filter(!(Cancer_type%in%bad_types$Cancer_type))
pd <- pd %>% select(Cancer_type,Plasma_volume_estimated) %>% distinct() %>% group_by(Cancer_type) %>% count() %>% filter(n>1) %>% select(-n) %>% inner_join(pd)

p <- ggplot(pd,aes(y=value,x=Plasma_volume_estimated)) 
p <- p + geom_boxplot(aes(fill=Plasma_volume_estimated),color="#525252",alpha=0.25,outlier.shape=16,outlier.alpha=0.5,outlier.size=1,width=0.5,linewidth=0.5)
p <- p + facet_wrap(~Cancer_type,nrow=1)
p <- p + scale_y_continuous(unique(pd$axisLabel),breaks=c(-1:3),labels=10**c(-1:3),expand = expansion(mult = c(0.15, 0.1)))
p <- p + scale_color_manual(values=c("#000000","#727272"))
p <- p + scale_fill_manual(values=c("#000000","#727272"))
p <- p + stat_compare_means(label = "p.format", method = "wilcox.test")
p <- p + theme_minimal() %+replace%
  theme_cowplot(12) %+replace%
  theme(
    plot.title = element_text(hjust=0,size=7,face="bold",margin=margin(0,0,4,0)),
    plot.subtitle = element_text(hjust=0,size=6,margin=margin(0,0,4,0),lineheight = 1),
    legend.position = "none",
    axis.title.x=element_blank(),
    axis.title.y = element_text(size=6,angle=90),
    axis.text.y = element_text(size=5,margin=margin(0,0,0,4)),
    axis.text.x = element_text(size=5.5,angle=0,hjust=0.5,margin=margin(2,0,0,0),vjust=1))

conc_boxes2 <- plot_boxes("All samples","log_DNA_conc",cancertypes,TRUE)


# --- Example of tumor fraction and lymphoma (Panel E) ---
statIn <- "TF"
# pdAll is already filtered to First_Sample
pd_tf <- pdAll %>% filter(stat=="TF" & !is.na(value))

# Build the groups manually as per original logic
pd_tf <- pd_tf %>% mutate(xset=Cancer_status) %>% filter(Cancer_status=="cancer")
pd_tf <- pd_tf %>% filter(species=="dog" & Cancer_type!="Lymphoma" & Cancer_status=="cancer") %>% mutate(xset="cancer (no_lymphomas)") %>% bind_rows(pd_tf)
pd_tf <- pd_tf %>% mutate(xset=paste(species,xset)) %>% distinct()
pd_tf <- pd_tf %>% select(xset,xpos) %>% group_by(xset) %>% summarize(xpos=min(xpos)) %>% inner_join(pd_tf %>% select(-xpos))

# Fix Labels
xlabels <- pd_tf %>% group_by(xset,xpos) %>% count() %>% mutate(xlabel=paste(str_replace(xset," \\(","\n("),"\nN=",n,sep="")) %>% select(-n)
xlabels <- xlabels %>% mutate(xlabel=str_replace_all(xlabel,"\\ ","\n"))
xlabels <- xlabels %>% mutate(xlabel=str_replace_all(xlabel,"\\_"," "))

min_max <- pd_tf %>% ungroup() %>% group_by(xset,xpos) %>% summarize(min=min(value),max=max(value))
min_max <- min_max %>% left_join(xlabels)
pd_tf <- pd_tf %>% mutate(xset=paste(Cancer_status,xset))
statstrings <- pd_tf %>% group_by(xset) %>% summarize(n=n(),mean=mean(value),sd=sd(value)) %>% mutate(statstr=paste(round(mean,3),"\n+/-",round(sd,3),sep=""))

xbreaks <- c("cancer dog cancer","cancer dog cancer (no_lymphomas)","cancer human cancer")
xlabels <- c("Dog cancers","Dog cancers\n(no lymphomas)","Human cancers")

statstrings <- statstrings %>% left_join(tibble(xset=xbreaks,labels=xlabels,order=c(1:3))) %>% mutate(labels=paste(labels,"\nN=",n,sep=""))
xlabels <- statstrings %>% arrange(order) %>% pull(labels)
compares <- list(c(xbreaks[1:2]),c(xbreaks[2:3])) 

# --- Calculate Stats for CSV Only (Panel E) ---
# Comparison 1
g1_name <- "cancer dog cancer"
g2_name <- "cancer dog cancer (no_lymphomas)"
d_sub <- pd_tf %>% filter(xset %in% c(g1_name, g2_name)) %>% mutate(xset = factor(xset))

p_wilcox_1 <- tryCatch(wilcox.test(value ~ xset, data=d_sub)$p.value, error=function(e) NA)
p_perm_t_1 <- tryCatch({
  .perm_p_indep(d_sub$value, d_sub$xset, B = 10000)
}, error=function(e) NA)

# Comparison 2
g3_name <- "cancer human cancer"
d_sub2 <- pd_tf %>% filter(xset %in% c(g2_name, g3_name)) %>% mutate(xset = factor(xset))

p_wilcox_2 <- tryCatch(wilcox.test(value ~ xset, data=d_sub2)$p.value, error=function(e) NA)
p_perm_t_2 <- tryCatch({
  .perm_p_indep(d_sub2$value, d_sub2$xset, B = 10000)
}, error=function(e) NA)

# Write to CSV with wilcox and perm p-values
stat.test.tf.csv <- tibble(
  group1 = c(xbreaks[1], xbreaks[2]),
  group2 = c(xbreaks[2], xbreaks[3]),
  wilcox_p = c(p_wilcox_1, p_wilcox_2),
  perm_t_p = c(p_perm_t_1, p_perm_t_2)
) %>%
  mutate(
    wilcox_p = format_p_csv(wilcox_p),
    perm_t_p = format_p_csv(perm_t_p)
  )

write_lines(c("", "Panel E Specific Comparisons (Wilcoxon & Permuted T-test)"), outputfile, append=TRUE)
write.table(stat.test.tf.csv, file=outputfile, append=TRUE, sep=",", row.names=FALSE)


# --- Plot Panel E (Wilcoxon) ---
rangey <- max(pd_tf$value)-min(pd_tf$value)
miny <- min(pd_tf$value) -  (rangey/20)
pd_tf$xset <- factor(pd_tf$xset,levels=xbreaks)

p <- ggplot(pd_tf,aes(y=value,x=xset)) 
p <- p + geom_text(aes(label=statstr),color="#cb181d",y=miny, data=statstrings,size=1.75,lineheight=0.9,hjust=0.5,vjust=1)
p <- p + geom_boxplot(fill="#cb181d",color="#525252",alpha=0.25,outlier.shape=16,outlier.alpha=0.5,outlier.size=1,width=0.5,linewidth=0.5)
p <- p + ggtitle(unique(pd_tf$statName),subtitle="First samples (Wilcoxon)")

p <- p + stat_compare_means(method="wilcox.test",comparisons=compares,size=2.5)

if (statIn=="log_DNA_conc"){
  p <- p + scale_y_continuous(unique(pd_tf$axisLabel),breaks=c(-1:3),labels=10**c(-1:3),expand = expansion(mult = c(0.15, 0.2)))
} else {
  p <- p + scale_y_continuous(unique(pd_tf$axisLabel),expand = expansion(mult = c(0.15, 0.2)))
}
p <- p + scale_x_discrete("",breaks=xbreaks,labels=xlabels) 
p <- p + scale_shape_manual(values=c(21,19))
p <- p + theme_minimal() %+replace%
  theme_cowplot(12) %+replace%
  theme(
    plot.title = element_text(hjust=0,size=7,face="bold",margin=margin(0,0,4,0)),
    plot.subtitle = element_text(hjust=0,size=6,margin=margin(0,0,4,0),lineheight = 1),
    legend.position = "none",
    axis.title.x=element_blank(),
    axis.title.y = element_text(size=6,angle=90),
    axis.text.y = element_text(size=5,margin=margin(0,0,0,4)),
    axis.text.x = element_text(size=6,angle=0,hjust=0.5,margin=margin(2,0,0,0),vjust=1))
p_TF_example <- p


# --- Assemble figure ---
row_table <- plot_grid(
  ft_plot,
  labels=c("A. Anova analysis on liquid biopsy metrics by species & cancer type (first samples)"),
  ncol=1,
  label_size = 12,
  hjust=0
)
row_boxplots1 <- plot_grid(tf_boxes,frag_boxes,labels=LETTERS[2:3],ncol=2,label_size = 12,rel_widths = c(1,1)) 
row_boxplots2 <- plot_grid(conc_boxes,p_TF_example,labels=LETTERS[4:5],ncol=2,label_size = 12,rel_widths = c(1.75,1)) 

pdfname <- paste(pdfdir,"Fig_S3.pdf",sep="")
grid <- plot_grid(row_table,row_boxplots1,row_boxplots2,labels = NULL,ncol=1) 
ggsave(grid,filename=pdfname,width=6.5,height=9)
print(paste("Made",pdfname))

# --- Sample metadata ---
sample_metadata_file <- paste0(pdfdir, "Fig_S3.sample_metadata.csv")
write.csv(pdAll, file = sample_metadata_file, row.names = FALSE)
message("Wrote sample metadata to: ", sample_metadata_file)