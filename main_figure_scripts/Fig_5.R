library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(lme4)
library(lmerTest)
library(patchwork)
if (!requireNamespace("robustbase", quietly = TRUE)) install.packages("robustbase")
library(robustbase)
# This package automatically finds the root folder of the downloaded github project.
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ==============================================================================
# 1. SETUP & FILE PATHS
# ==============================================================================

# Define the input file (pointing to the 'data' folder in the repo)
infile <- here("data", "BDS.txt")

# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

# Mapping of metric codes to remove spaces
rename_stats <- tibble(
  order     = c(2, 1, 3),
  statistic = c("Tumor fraction","log10 cfDNA concentration","Fragment size ratio"),
  stat      = c("TF","log_DNA_conc","Fragment_size_ratio")
) %>% arrange(order)

# --- Helper: Formatter for CSV ---
format_p_csv <- function(x) {
  sapply(x, function(val) {
    if (is.na(val) || is.infinite(val)) return("NA")
    if (val == 0) return("< 2.2e-16") 
    if (val < 0.001) {
      formatC(val, format = "e", digits = 4)
    } else {
      formatC(val, format = "f", digits = 4)
    }
  })
}

## ------------------------------------------------------------------
## Permutation helper: Supports Mixed Models
## ------------------------------------------------------------------

nperm_anova <- 100      # number of permutations
trim_frac   <- 0.05     # fraction trimmed

perm_for_model <- function(dsubset, res, formula, nperm = 100, trim_frac = 0.05) {
  out_list <- list()
  stats <- unique(dsubset$stat)
  
  # Check if formula implies a mixed model
  is_mixed <- grepl("\\|", as.character(formula)[3])
  
  for (s in stats) {
    df <- dsubset %>% filter(stat == s)
    if (nrow(df) < 5) next
    
    # Effects we actually care about
    effs <- res %>% filter(stat == s) %>% pull(Effect)
    if (length(effs) == 0) next
    
    # --- Get Observed F ---
    an_tbl <- tryCatch({
      if (is_mixed) {
        m_obs <- suppressMessages(lmer(formula, data = df))
        at <- suppressMessages(anova(m_obs))
        at <- as.data.frame(at) %>% rownames_to_column("Effect") %>% as_tibble()
        at %>% rename(F = `F value`)
      } else {
        anova_test(data = df, formula = formula) %>% as_tibble()
      }
    }, error = function(e) { NULL })
    
    if(is.null(an_tbl)) next
    
    F_obs  <- an_tbl$F[match(effs, an_tbl$Effect)]
    k      <- length(effs)
    ge_full <- numeric(k)
    ge_trim <- numeric(k)
    
    # --- Determine trimmed subset ---
    n <- nrow(df)
    keep_trim <- integer(0)
    if (n >= 5) {
      q_low  <- quantile(df$value, trim_frac, na.rm = TRUE)
      q_high <- quantile(df$value, 1 - trim_frac, na.rm = TRUE)
      keep_trim <- which(df$value >= q_low & df$value <= q_high)
      
      min_obs <- if(is_mixed) length(fixef(m_obs)) + 2 else length(coef(lm(formula, data=df))) + 1
      if (length(keep_trim) <= min_obs) {
        keep_trim <- integer(0)
      }
    }
    
    if (k > 0) {
      for (i in seq_len(nperm)) {
        df_perm <- df
        df_perm$value <- sample(df$value) # Simple shuffle of Y
        
        # --- Run Permuted Model ---
        F_perm_full <- rep(NA, k)
        
        if (is_mixed) {
          tryCatch({
            m_perm <- suppressMessages(lmer(formula, data = df_perm))
            if(!is.null(m_perm)) {
              an_perm <- suppressMessages(anova(m_perm))
              an_perm <- as.data.frame(an_perm) %>% rownames_to_column("Effect")
              F_perm_full <- an_perm$`F value`[match(effs, an_perm$Effect)]
            }
          }, error = function(e) { NULL })
        } else {
          tryCatch({
            an_perm_full <- anova_test(data = df_perm, formula = formula) %>% as_tibble()
            F_perm_full  <- an_perm_full$F[match(effs, an_perm_full$Effect)]
          }, error = function(e) { NULL })
        }
        
        ge_full <- ge_full + as.numeric(!is.na(F_perm_full) & !is.na(F_obs) & F_perm_full >= F_obs)
        
        # --- Run Trimmed Permuted Model ---
        if (length(keep_trim) > 0) {
          df_perm_trim <- df_perm[keep_trim, , drop = FALSE]
          F_perm_trim <- rep(NA, k)
          
          if (is_mixed) {
            tryCatch({
              m_perm_trim <- suppressMessages(lmer(formula, data = df_perm_trim))
              if(!is.null(m_perm_trim)) {
                an_perm_trim <- suppressMessages(anova(m_perm_trim))
                an_perm_trim <- as.data.frame(an_perm_trim) %>% rownames_to_column("Effect")
                F_perm_trim <- an_perm_trim$`F value`[match(effs, an_perm_trim$Effect)]
              }
            }, error = function(e) { NULL })
          } else {
            tryCatch({
              an_perm_trim <- anova_test(data = df_perm_trim, formula = formula) %>% as_tibble()
              F_perm_trim  <- an_perm_trim$F[match(effs, an_perm_trim$Effect)]
            }, error = function(e) { NULL })
          }
          
          ge_trim <- ge_trim + as.numeric(!is.na(F_perm_trim) & !is.na(F_obs) & F_perm_trim >= F_obs)
        }
      }
    }
    
    p_perm      <- if (k > 0) ge_full / nperm else numeric(0)
    p_perm_trim <- if (length(keep_trim) > 0 && k > 0) ge_trim / nperm else rep(NA_real_, k)
    
    out_list[[length(out_list) + 1]] <- tibble(
      stat        = s,
      Effect      = effs,
      p_perm      = p_perm,
      p_perm_trim = p_perm_trim
    )
  }
  
  if (length(out_list) == 0) {
    tibble(stat = character(0), Effect = character(0), p_perm = double(0), p_perm_trim = double(0))
  } else {
    bind_rows(out_list)
  }
}

## ------------------------------------------------------------------
## Load and merge data
## ------------------------------------------------------------------

if(file.exists(infile)) {
  # 1. LOAD DATA (UPDATED: Single long-form spreadsheet)
  raw <- as_tibble(read.delim(
    infile,
    header = TRUE,
    check.names = FALSE, # Preserve T_simple / T_circulation
    na.strings = c("NA","N/A","","UNK","?")
  ))
  
  # 2. CLEAN & PROCESS
  # Convert types, handle Log DNA, Rename FSR
  processed_data <- raw %>%
    mutate(
      Dog_ID = as.factor(Dog_ID),
      Sample_ID = trimws(Sample_ID),
      Blood_draw_site = trimws(Blood_draw_site),
      # Calculate Log DNA Conc
      log_DNA_conc = log10(DNA_conc),
      # Ensure numeric distances
      T_simple = as.numeric(T_simple),
      T_circulation = as.numeric(T_circulation)
    ) %>%
    # Rename FSR to Fragment_size_ratio
    rename(Fragment_size_ratio = FSR)
  
  # 3. PREPARE ANALYSIS DATASET (Wide Format for Regression)
  # One row per sample, with columns for metrics and BOTH distance types
  # Pivot metrics longer so we can iterate over them (TF, log_DNA_conc, FSR)
  d_expanded <- processed_data %>%
    select(
      Sample_ID, Dog_ID, Blood_draw_site,
      T_simple, T_circulation,
      TF, log_DNA_conc, Fragment_size_ratio
    ) %>%
    pivot_longer(
      cols = c(TF, log_DNA_conc, Fragment_size_ratio),
      names_to = "stat",
      values_to = "value"
    ) %>%
    filter(!is.na(value))
  
  # 4. PREPARE PLOT DATA (Long Format for Boxplots)
  # Boxplots only need metrics and site, d_expanded works for this too
  # (pSite_conc, etc. will filter by 'stat')
  
  cat("--------------------------------------------------------\n")
  cat("DIAGNOSTICS:\n")
  cat("Number of unique samples:", n_distinct(d_expanded$Sample_ID), "\n")
  cat("Number of rows in Analysis Dataset:", nrow(d_expanded), "\n")
  cat("--------------------------------------------------------\n")
  
  ## Export expanded table (User verification)
  analysis_table_file <- paste(pdfdir, "Table.Fig_5.analysis_data.csv", sep = "")
  write.table(d_expanded, file = analysis_table_file, sep = ",", row.names = FALSE)
  
} else {
  stop(paste("Input file not found:", infile))
}

## ------------------------------------------------------------------
## Plot themes
## ------------------------------------------------------------------

theme_dist <- function(){
  theme_cowplot(12) %+replace%
    theme(
      plot.title    = element_text(hjust=0,size=9,face="bold",
                                   margin=margin(0,0,4,0)),
      plot.subtitle= element_text(hjust=0,size=9,
                                  margin=margin(0,0,4,0),lineheight = 1),
      axis.title    = element_text(size=9),
      axis.text     = element_text(size=8),
      panel.spacing= unit(1.2, "lines"),
      legend.position="NONE"
    )
}

theme_site <- function(){
  theme_cowplot(12) %+replace%
    theme(
      plot.title    = element_text(hjust=0,size=9,face="bold",
                                   margin=margin(0,0,4,0)),
      plot.subtitle= element_text(hjust=0,size=9,
                                  margin=margin(0,0,4,0),lineheight = 1),
      axis.title    = element_text(size=9),
      axis.text     = element_text(size=8),
      axis.ticks.x = element_blank(),
      panel.spacing= unit(1.2, "lines"),
      legend.justification = "center",
      legend.text = element_text(size=8),
      legend.title = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.margin = margin(0,0,0,0)
    )
}

## ------------------------------------------------------------------
## Plots
## ------------------------------------------------------------------

my_comparisons <- list(
  c("Jugular","Cephalic"),
  c("Jugular","Saphenous"),
  c("Cephalic","Saphenous")
)

# UPDATED: Dist vs site plot
# We need to pivot T_simple and T_circulation to long format for this plot
pdDist <- processed_data %>%
  select(Sample_ID, Blood_draw_site, T_simple, T_circulation, Dog_ID) %>%
  pivot_longer(cols = c(T_simple, T_circulation), 
               names_to = "distance_type", 
               values_to = "distance") %>%
  filter(!is.na(distance)) %>%
  mutate(distance_type = dplyr::recode(distance_type, 
                                       "T_simple" = "Simple distance", 
                                       "T_circulation" = "Circulation distance"))

pdDist$Blood_draw_site <- factor(
  pdDist$Blood_draw_site,
  levels=c("Jugular","Cephalic","Saphenous")
)

# Calculate distinct dogs per site for distance plot
dist_counts <- pdDist %>%
  group_by(Blood_draw_site) %>%
  summarise(n = n_distinct(Dog_ID), .groups = 'drop') %>%
  mutate(label = case_when(
    Blood_draw_site == "Saphenous" ~ paste0("Saph-\nenous\n(n=", n, ")"),
    TRUE ~ paste0(Blood_draw_site, "\n(n=", n, ")")
  ))

dist_labels <- setNames(dist_counts$label, dist_counts$Blood_draw_site)

p <- ggplot(pdDist,aes(x=Blood_draw_site, y=distance))
p <- p + geom_boxplot(aes(fill=Blood_draw_site),
                      outlier.shape = NA,width=0.5,alpha=0.1,linewidth=0.5)
p <- p + geom_jitter(aes(color=Blood_draw_site,shape=Blood_draw_site),
                     size=1.5,alpha=0.75,height=0,width=0.1)
p <- p + stat_compare_means(comparisons = my_comparisons,size=2.5)
p <- p + facet_wrap(~distance_type, scales="free_y")
p <- p + scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))
p <- p + scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))
p <- p + scale_shape_manual(values=c(16:18))
p <- p + scale_x_discrete("blood draw site",
                          breaks=names(dist_labels),
                          labels=dist_labels)
p <- p + ggtitle("Distances by vein sampled")
p <- p + theme_site()
p <- p + theme(legend.position = "none")
p_dist_v_site <- p

# Metrics by Site Plot
plot_corrbox_plot <- function(data, nameIn, titleIn){
  subset <- data %>%
    filter(stat == nameIn) %>%
    left_join(rename_stats, by = "stat") %>%
    select(Dog_ID, Sample_ID, Blood_draw_site, value, stat, statistic) %>%
    distinct()
  
  subset$Blood_draw_site <- factor(subset$Blood_draw_site, levels=c("Jugular","Cephalic","Saphenous"))
  
  # Calculate distinct dogs per site for dynamic labeling
  n_counts <- subset %>%
    group_by(Blood_draw_site) %>%
    summarise(n = n_distinct(Dog_ID), .groups = 'drop') %>%
    mutate(label = case_when(
      Blood_draw_site == "Saphenous" ~ paste0("Saph-\nenous\n(n=", n, ")"),
      TRUE ~ paste0(Blood_draw_site, "\n(n=", n, ")")
    ))
  
  site_labels <- setNames(n_counts$label, n_counts$Blood_draw_site)
  
  p <- ggplot(subset,aes(x=Blood_draw_site,y=value))
  p <- p + geom_boxplot(aes(fill=Blood_draw_site),
                        outlier.shape = NA,width=0.5,alpha=0.1,linewidth=0.5)
  p <- p + geom_line(aes(group=Dog_ID),color="#525252",linewidth=0.25)
  p <- p + geom_point(aes(color=Blood_draw_site,shape=Blood_draw_site),
                      size=1.5,alpha=0.8)
  p <- p + scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))
  p <- p + scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))
  p <- p + scale_shape_manual(values=c(16:18))
  
  p <- p + scale_y_continuous(
    unique(subset$statistic),
    expand = expansion(mult = c(0.25, 0.1))
  )
  
  p <- p + scale_x_discrete("blood draw site",
                            breaks=names(site_labels),
                            labels=site_labels)
  p <- p + ggtitle(titleIn)
  p <- p + theme_site()
  p <- p + theme(axis.title.x = element_blank())
  p
}

## ------------------------------------------------------------------
## ANOVA + permutation checks (UPDATED MODEL FORMULA)
## ------------------------------------------------------------------

# Output stats file for this figure
outputfile <- paste(pdfdir,"Table.Fig_5.stats.csv",sep="")
write_lines(
  "Statistical results for analyses of blood draw site and separated distances",
  outputfile,
  append = FALSE
)

# Function to run model suite on a dataset
run_site_model_suite <- function(dataset, set_label, set_id) {
  # Note: dataset already pivoted long by metric (stat, value) 
  # but has T_simple and T_circulation as columns
  dsubset <- dataset %>% 
    filter(!is.na(value), !is.na(T_simple), !is.na(T_circulation), 
           !is.na(Blood_draw_site), !is.na(Dog_ID))
  
  ntotal <- dsubset %>% group_by(stat) %>% count()
  nsamples <- length(unique(dsubset$Sample_ID))
  ndogs    <- length(unique(dsubset$Dog_ID))
  sampleNumbers <- paste(ndogs," dogs",sep="")
  
  res_list <- list()
  stats_list <- unique(dsubset$stat)
  
  for(s in stats_list) {
    df_s <- dsubset %>% filter(stat == s)
    
    # 1. Main Mixed Model (Using T_simple and T_circulation directly)
    m <- lmer(value ~ Blood_draw_site + T_simple + T_circulation + (1|Dog_ID), data = df_s)
    
    an_tbl <- anova(m) %>% as.data.frame() %>% rownames_to_column("Effect") %>% as_tibble()
    an_tbl <- an_tbl %>% 
      rename(F = `F value`, p = `Pr(>F)`) %>%
      mutate(
        stat = s, 
        DFn = NumDF, 
        DFd = DenDF, 
        ges = (F * DFn) / (F * DFn + DFd)
      )
    
    # 2. Robust Approximation (lmrob on Fixed Effects)
    an_tbl$p_robust <- NA
    f_robust <- as.formula("value ~ Blood_draw_site + T_simple + T_circulation")
    
    try({
      rob_mod <- lmrob(f_robust, data = df_s, setting="KS2014")
      coefs <- summary(rob_mod)$coefficients
      
      for(i in 1:nrow(an_tbl)) {
        eff <- an_tbl$Effect[i]
        matches <- rownames(coefs)[grep(paste0("^", eff), rownames(coefs))]
        if(length(matches) > 0) {
          an_tbl$p_robust[i] <- min(coefs[matches, "Pr(>|t|)"], na.rm=TRUE)
        }
      }
    }, silent=TRUE)
    
    res_list[[s]] <- an_tbl
  }
  res <- bind_rows(res_list)
  
  perm <- perm_for_model(
    dsubset,
    res,
    value ~ Blood_draw_site + T_simple + T_circulation + (1|Dog_ID),
    nperm = nperm_anova,
    trim_frac = trim_frac
  )
  res <- res %>% left_join(perm, by = c("stat","Effect"))
  
  out <- res %>%
    mutate(
      setIn         = set_id,
      modelName     = "Blood_draw_site+T_simple+T_circ+(1|Dog_ID)",
      sampleNumbers = sampleNumbers,
      set           = paste(set_label, sampleNumbers, sep="; ")
    ) %>%
    left_join(ntotal, by = "stat")
  
  return(out)
}

# --- Analysis 1: All Dogs (Full Cohort) ---
res1 <- run_site_model_suite(d_expanded, "All samples", 1)

# --- Analysis 2: Shedders Only (Dogs with at least one TF > 0.03) ---
dogs_no_shed <- d_expanded %>% 
  filter(stat == "TF") %>%
  group_by(Dog_ID) %>%
  summarize(max_tf = max(value, na.rm=TRUE)) %>%
  filter(max_tf <= 0.03) %>%
  pull(Dog_ID)

d_shedders <- d_expanded %>% filter(!Dog_ID %in% dogs_no_shed)
res2 <- run_site_model_suite(d_shedders, "Dogs with TF > 0.03", 2)

# Combine Results
anovaOut <- bind_rows(res1, res2)

rename_effects <- tibble(
  Effect = c("Blood_draw_site","T_simple","T_circulation"),
  label  = c("Blood draw site","Simple distance","Circulation distance")
) %>% mutate(order = row_number())

anovaOut <- anovaOut %>%
  left_join(rename_stats, by = "stat") %>%
  mutate(
    p.adj = p.adjust(p, method = "fdr"),
    sig   = if_else(p <= 0.05, "significant", "not significant")
  )

anovaOut$statistic <- factor(
  anovaOut$statistic,
  levels = rename_stats %>% arrange(order) %>% pull(statistic)
)

rename_effects <- rename_effects %>% arrange(order)
anovaOut <- anovaOut %>% filter(Effect %in% rename_effects$Effect)
anovaOut$Effect <- factor(anovaOut$Effect, levels = rename_effects$Effect)

anovaOut <- anovaOut %>%
  mutate(
    pstr = if_else(
      p < 0.001,
      paste("p=", format(p, digits=3, scientific=TRUE), sep=""),
      paste("p=", as.character(signif(p, 3)), sep="")
    )
  )

anovaOut <- anovaOut %>%
  mutate(
    flag_outlier_driven = !is.na(p.adj)  & (p.adj <= 0.05) &
      !is.na(p_perm) & (p_perm <= 0.05) &
      !is.na(p_perm_trim) & (p_perm_trim > 0.05),
    flag_perm = case_when(
      !is.na(p.adj) & (p.adj <= 0.05) & !is.na(p_perm) & (p_perm > 0.05)  ~ TRUE,
      !is.na(p.adj) & (p.adj > 0.05)  & !is.na(p_perm) & (p_perm <= 0.05) ~ TRUE,
      TRUE ~ FALSE
    )
  )

# Write full CSV
csv_output <- anovaOut %>%
  select(set, modelName, sampleNumbers, statistic, Effect, DFn, DFd, F, ges, 
         p, p.adj, p_perm, p_robust, n) %>%
  mutate(across(starts_with("p"), format_p_csv))

write.table(csv_output, file=outputfile, append=TRUE, sep=",", row.names=FALSE)

## ------------------------------------------------------------------
## REVISED ANOVA PLOT + LAYOUT (DUAL PDF OUTPUT)
## ------------------------------------------------------------------

theme_anova_panel <- function(){
  theme_cowplot(12) %+replace%
    theme(
      axis.line.x  = element_line(linewidth=0.25),
      axis.ticks.x = element_line(linewidth=0.25),
      axis.line.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title    = element_text(size=9), 
      axis.text.x   = element_text(size=8),  
      axis.text.y   = element_text(size=8,hjust=1,vjust=0.5, margin = margin(r = 2)),
      panel.border = element_rect(color = "#525252", fill = NA, linewidth = 0.5),
      panel.grid.major.x = element_line(color="#bdbdbd",linewidth = 0.15),
      legend.position = "none"
    )
}

create_single_anova <- function(data_anova, stat_code, show_y_labels = FALSE) {
  plot_data <- data_anova %>% filter(stat == stat_code)
  
  p <- ggplot(plot_data, aes(x=ges, y=Effect)) +
    geom_bar(aes(fill=sig), stat="identity", width=0.5, show.legend = FALSE) +
    geom_vline(xintercept = 0, color="#969696", linewidth=0.2) +
    geom_text(aes(label=pstr, color=sig),
              hjust=0, size=2.5, nudge_x = 0.02, show.legend = FALSE) +
    scale_y_discrete("",
                     breaks=rename_effects$Effect,
                     labels=rename_effects$label) +
    scale_color_manual(values=c("#525252","#000000")) +
    scale_fill_manual(values=c("#969696","#000000")) +
    scale_x_continuous(
      "Anova effect size (ges)",
      limits = c(0, 1),
      expand = c(0, 0),
      breaks = c(0, 0.25, 0.50, 0.75, 1),
      labels = c("0", "", "0.50", "", "1")
    ) +
    coord_cartesian(clip = "off") +
    theme_anova_panel()
  
  if(!show_y_labels) {
    p <- p + theme(axis.text.y = element_blank())
  }
  
  return(p)
}

# --- MASTER FUNCTION TO GENERATE AND SAVE FIGURE ---
generate_and_save_figure <- function(metrics_data, anova_data, filename_suffix) {
  
  # 1. Generate Boxplots
  pSite_conc <- plot_corrbox_plot(metrics_data, "log_DNA_conc", "log10 cfDNA concentration\nby vein sampled")
  pSite_tf   <- plot_corrbox_plot(metrics_data, "TF", "Tumor fraction\nby vein sampled")
  pSite_frag <- plot_corrbox_plot(metrics_data, "Fragment_size_ratio", "Fragment size ratio\nby vein sampled")
  
  # 2. Generate ANOVA Plots
  anova_conc <- create_single_anova(anova_data, "log_DNA_conc", show_y_labels = TRUE)
  anova_tf   <- create_single_anova(anova_data, "TF", show_y_labels = TRUE)
  anova_frag <- create_single_anova(anova_data, "Fragment_size_ratio", show_y_labels = TRUE)
  
  # 3. Assemble
  col1 <- pSite_conc / anova_conc + plot_layout(heights = c(6, 1))
  col2 <- pSite_tf / anova_tf + plot_layout(heights = c(6, 1))
  col3 <- pSite_frag / anova_frag + plot_layout(heights = c(6, 1))
  
  fig_final <- (col1 | col2 | col3) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  # 4. Save
  fname <- paste0(pdfdir, "Fig_5", filename_suffix, ".pdf")
  ggsave(fig_final, filename=fname, width=10, height=5)
  print(paste("Saved:", fname))
}

# --- EXECUTE FOR BOTH VERSIONS ---

# Version 1: All Samples
generate_and_save_figure(
  metrics_data = d_expanded, 
  anova_data = anovaOut %>% filter(setIn == 1), 
  filename_suffix = "All"
)

# Version 2: Filtered (TF > 0.03)
generate_and_save_figure(
  metrics_data = d_expanded %>% filter(!Dog_ID %in% dogs_no_shed), 
  anova_data = anovaOut %>% filter(setIn == 2), 
  filename_suffix = "Filtered"
)

# Save the Distance vs Site plot separately (Using Long Data for Faceted Plot)
ggsave(p_dist_v_site, filename=paste(pdfdir,"Fig_S15_dist.pdf",sep=""), width=6, height=4)