library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(lubridate)
library(coin) 
library(lme4)     # For Mixed Models (Longitudinal)
library(lmerTest) # For p-values in LMM
library(patchwork)
library(robustbase) # For robust regression (lmrob)
library(WRS2)       # For robust paired tests (yuend)
# This package automatically finds the root folder of the downloaded github project.
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ==============================================================================
# 1. SETUP & FILE PATHS
# ==============================================================================

# Define the input file
infile <- here("data", "BB_Supplementary_Data_1.txt")

# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

outputfile <- paste0(pdfdir, "Table.Fig_S17.stats.csv")
write_lines(c("Statistical results for analyses presented in Fig. S17"), outputfile, append=FALSE)

rename_stats <- tibble(
  order=c(2,1,3),
  statName=c("Tumor fraction","cfDNA concentration","Fragment size ratio"),
  axisLabel=c("Tumor fraction","cfDNA concentration (ng/mL)","Fragment size ratio"),
  stat=c("TF","log_DNA_conc","Fragment_size_ratio") 
) %>% arrange(order)

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

format_p_plot <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ formatC(p, format = "e", digits = 1),
    TRUE ~ formatC(p, format = "f", digits = 3)
  )
}

# --- Data Loading ---
if(file.exists(infile)) {
  raw <- as_tibble(read.csv(infile, header=T, sep="\t", na.strings = c("NA","N/A","","UNK","?"))) 
  raw <- raw %>% mutate(Sample_type=str_remove(Sample_type,"_plasma"))
  raw <- raw %>% mutate(Date_of_sample=str_replace(Date_of_sample," \\(approx\\)",""))
  
  # Parse Times
  raw <- raw %>% filter(Cohort=="Time_of_day")
  
  # Parse Times
  raw <- raw %>% 
    filter(!is.na(Date_of_sample) & !is.na(Blood_draw_time)) %>% 
    mutate(Blood_draw_time=paste(Date_of_sample, Blood_draw_time)) %>%
    mutate(time=parse_date_time(Blood_draw_time, "m/d/y HMS Op", tz = "America/New_York")) %>%
    mutate(minutes=hour(time)*60 + minute(time))
  
  # Filter Duplicates
  raw <- raw %>% 
    group_by(ID, Date_of_sample, time) %>% 
    distinct() %>% 
    group_by(ID, Date_of_sample) %>% 
    mutate(n_samples = n()) %>%
    filter(n_samples > 1) %>% 
    ungroup() %>%
    mutate(log_DNA_conc = log10(as.numeric(DNA_conc)),
           TF = as.numeric(TF),
           Fragment_size_ratio = as.numeric(Fragment_size_ratio))
  
  # Filter for sufficient time range (>60 min)
  raw <- raw %>% 
    group_by(ID, Date_of_sample) %>% 
    summarize(diff = max(minutes) - min(minutes), .groups="drop") %>% 
    inner_join(raw, by=c("ID", "Date_of_sample")) %>%
    mutate(Cancer_status = if_else(Cancer_type != "Healthy", "cancer", "healthy"))
  
  # Prepare Long Data
  pdAll_base <- raw %>% 
    select(ID, Cancer_status, Cancer_type, log_DNA_conc, TF, Fragment_size_ratio, time, minutes) %>%
    pivot_longer(c(log_DNA_conc, TF, Fragment_size_ratio), names_to="stat", values_to="value") %>% 
    filter(!is.na(value)) %>%
    group_by(ID, stat) %>% 
    mutate(maxtime = max(minutes), mintime = min(minutes), lapse = minutes - mintime) %>%
    ungroup() %>%
    left_join(rename_stats, by="stat")
  
  # Normalize to Start Value (for relative plots)
  mins <- pdAll_base %>% 
    group_by(ID, stat) %>% 
    filter(minutes == mintime) %>%
    summarize(start_value = mean(value), .groups="drop") 
  
  pdAll_base <- pdAll_base %>% left_join(mins, by=c("ID", "stat")) %>%
    mutate(rel = if_else(stat == "log_DNA_conc", 10^value / 10^start_value, value / start_value)) %>%
    filter(start_value > 0)
  
} else {
  stop("Input file not found.")
}

# ==============================================================================
# ANALYSIS FUNCTION (Longitudinal + Paired)
# ==============================================================================

run_stime_analysis <- function(dataset, set_label) {
  
  # --- PART 1: LONGITUDINAL TREND (LMM) ---
  anova_stats <- dataset %>%
    group_by(stat) %>%
    group_modify(~ {
      d_sub <- .x
      
      # 1. Linear Mixed Model (Primary)
      # Controls for dog ID as random effect
      res_lmm <- tryCatch({
        m <- lmer(value ~ minutes + (1|ID), data = d_sub)
        coefs <- summary(m)$coefficients
        coefs["minutes", "Pr(>|t|)"]
      }, error=function(e) NA)
      
      # 2. Permutation Test (Robust to distribution)
      res_perm <- tryCatch({
        d_perm <- d_sub %>% mutate(ID = factor(ID)) %>% as.data.frame()
        as.numeric(coin::pvalue(coin::independence_test(value ~ minutes | ID, data = d_perm, distribution="asymptotic")))
      }, error=function(e) NA)
      
      # 3. Robust Linear Regression (Robust Check)
      # Uses MM-estimation to down-weight outliers (robustbase::lmrob)
      # Note: This ignores random effect of ID, serving as a robust population trend check
      res_robust <- tryCatch({
        rob_fit <- lmrob(value ~ minutes, data = d_sub, setting="KS2014")
        summary(rob_fit)$coefficients["minutes", "Pr(>|t|)"]
      }, error=function(e) NA)
      
      tibble(
        p_lmm = res_lmm,
        p_perm = res_perm,
        p_robust = res_robust,
        n_samples = nrow(d_sub),
        n_dogs = n_distinct(d_sub$ID)
      )
    }) %>%
    mutate(set = set_label)
  
  # Annotate dataset with stats for plotting
  dataset_labeled <- dataset %>% 
    left_join(anova_stats, by="stat") %>%
    mutate(plot_subtitle = paste("LMM p =", format_p_plot(p_lmm)))
  
  # --- PART 2: MORNING VS AFTERNOON (PAIRED) ---
  pdBox <- dataset %>% 
    select(ID, stat, time, value, maxtime, mintime, minutes) %>% 
    mutate(timeset = if_else(hour(time) < 12, "morning", "afternoon")) %>%
    # Keep only earliest morning and latest afternoon
    filter((timeset=="morning" & mintime==minutes) | (timeset=="afternoon" & maxtime==minutes)) %>%
    distinct(ID, stat, timeset, value) %>%
    group_by(ID, stat) %>% 
    filter(n() == 2) %>% 
    ungroup()
  
  pdBox$timeset <- factor(pdBox$timeset, levels=c("morning", "afternoon"))
  
  paired_stats <- pdBox %>%
    group_by(stat) %>%
    group_modify(~ {
      d_wide <- .x %>% pivot_wider(names_from=timeset, values_from=value)
      
      # Parametric
      p_ttest <- tryCatch(t.test(d_wide$morning, d_wide$afternoon, paired=TRUE)$p.value, error=function(e) NA)
      
      # Non-parametric
      p_wilcox <- tryCatch(wilcox.test(d_wide$morning, d_wide$afternoon, paired=TRUE)$p.value, error=function(e) NA)
      
      # Permutation
      d_long <- .x %>% mutate(ID = factor(ID), timeset = factor(timeset)) %>% as.data.frame()
      p_perm <- tryCatch({
        as.numeric(coin::pvalue(coin::symmetry_test(value ~ timeset | ID, data=d_long, distribution="asymptotic")))
      }, error=function(e) NA)
      
      # Robust (Trimmed Means, Yuen's Test)
      p_robust <- tryCatch({
        WRS2::yuend(d_wide$morning, d_wide$afternoon, tr=0.2)$p.value
      }, error=function(e) NA)
      
      tibble(
        p_paired_t = p_ttest,
        p_paired_wilcox = p_wilcox,
        p_perm = p_perm,
        p_robust_trimmed = p_robust,
        n_pairs = nrow(d_wide)
      )
    }) %>%
    mutate(set = set_label)
  
  # Annotate box dataset
  pdBox_labeled <- pdBox %>% 
    left_join(paired_stats, by="stat") %>%
    mutate(subtitle_box = paste("Paired Wilcoxon p =", format_p_plot(p_paired_wilcox)))
  
  return(list(
    long_stats = anova_stats,
    paired_stats = paired_stats,
    data_long = dataset_labeled,
    data_box = pdBox_labeled
  ))
}

# ==============================================================================
# EXECUTE ANALYSES
# ==============================================================================

# 1. Set 1: All Dogs
res1 <- run_stime_analysis(pdAll_base, "All dogs")

# 2. Set 2: Filtered (TF > 0.03)
# Identify dogs where ALL TF measurements are <= 0.03
dogs_no_shed <- pdAll_base %>% 
  filter(stat == "TF") %>%
  group_by(ID) %>%
  summarize(max_tf = max(value, na.rm=TRUE)) %>%
  filter(max_tf <= 0.03) %>%
  pull(ID)

pdAll_filtered <- pdAll_base %>% filter(!ID %in% dogs_no_shed)
res2 <- run_stime_analysis(pdAll_filtered, "Dogs with TF > 0.03")

# ==============================================================================
# SAVE STATS TO CSV
# ==============================================================================

# Combine Longitudinal Stats
long_combined <- bind_rows(res1$long_stats, res2$long_stats) %>%
  select(set, stat, n_dogs, n_samples, p_lmm, p_perm, p_robust) %>%
  mutate(across(starts_with("p"), format_p_csv))

write_lines(c("", "Part 1: Longitudinal Effect of Time (LMM vs Permutation vs Robust Regression)"), outputfile, append=TRUE)
write.table(long_combined, file=outputfile, append=TRUE, sep=",", row.names=FALSE)

# Combine Paired Stats
paired_combined <- bind_rows(res1$paired_stats, res2$paired_stats) %>%
  select(set, stat, n_pairs, p_paired_wilcox, p_perm, p_robust_trimmed, p_paired_t) %>%
  mutate(across(starts_with("p"), format_p_csv))

write_lines(c("", "Part 2: Morning vs Afternoon (Paired Wilcoxon vs Permutation vs Robust Yuen vs T-test)"), outputfile, append=TRUE)
write.table(paired_combined, file=outputfile, append=TRUE, sep=",", row.names=FALSE)

# ==============================================================================
# PLOTTING
# ==============================================================================

theme_replicates <- function(){ 
  theme_cowplot(12) %+replace%
    theme(
      plot.title = element_text(hjust=0,size=9,face="bold",margin=margin(0,0,4,0)),
      plot.subtitle = element_text(hjust=0,size=8,margin=margin(0,0,4,0),lineheight = 1),
      axis.title = element_text(size=9),
      axis.text = element_text(size=8,margin=margin(2,0,0,4)),
      legend.position="none"
    )
}

# --- Plot Generator Function ---
generate_stime_figure <- function(data_long, data_box, filename_suffix) {
  
  # Setup Axis
  xmin <- floor(min(data_long$minutes)/60)/2
  xmax <- floor(max(data_long$minutes)/60)/2
  xaxis <- tibble(breaks=c(xmin:xmax)*2*60) %>%
    mutate(labels=breaks/60) %>%
    mutate(labels=if_else(labels<12, paste(labels,"am"), if_else(labels > 12, paste(labels-12,"pm"), paste(labels,"pm"))))
  
  # Colors/Shapes
  unique_ids <- unique(data_long$ID)
  colors <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3") 
  if(length(unique_ids) > 4) colors <- rep(colors, ceiling(length(unique_ids)/4))
  shapes <- c(15:17)
  shapes <- rep(shapes, ceiling(length(unique_ids)/3))
  
  # 1. Longitudinal Plots
  plot_dist_plot <- function(statIn){
    pd <- data_long %>% filter(stat==statIn) 
    statNameIn <- unique(pd$statName)
    subtitle <- unique(pd$plot_subtitle)
    
    p <- ggplot(pd, aes(x=minutes, y=rel)) +
      geom_point(aes(color=ID, shape=ID), alpha=0.75, size=1.5) +
      geom_line(aes(group=ID, color=ID), alpha=0.75) +
      ggtitle(statNameIn, subtitle=subtitle) +
      scale_x_continuous("time", breaks=xaxis$breaks, labels=xaxis$labels) +
      scale_color_manual(values=colors) +
      scale_shape_manual(values=shapes) +
      scale_y_continuous(paste("relative", statNameIn)) +
      theme_replicates() +
      theme(legend.position=c(0.02,0.98),
            legend.justification = c("left", "top"),
            legend.margin = margin(2,2,2,2),
            legend.text = element_text(size=6, vjust=0.5), # Increased font
            legend.key.size = unit(0.3, 'cm'),
            legend.title = element_blank(),
            legend.box.background = element_rect(color="#d9d9d9", size=0.5))
    return(p)
  }
  
  pDist_conc <- plot_dist_plot("log_DNA_conc")
  pDist_tf <- plot_dist_plot("TF")
  pDist_frag <- plot_dist_plot("Fragment_size_ratio")
  
  # 2. Paired Boxplots
  plot_timesets <- function(nameIn){
    subset <- data_box %>% filter(stat==nameIn) %>% distinct()
    statName <- rename_stats %>% filter(stat==nameIn) %>% pull(statName)
    subtitle <- unique(subset$subtitle_box)
    
    p <- ggplot(subset, aes(x=timeset, y=value)) +
      geom_boxplot(outlier.shape = NA, width=0.5, alpha=0.1, linewidth=0.5) +
      geom_line(aes(group=ID, color=ID), linewidth=0.25) +
      geom_point(aes(color=ID, shape=ID), size=1.5, alpha=0.75) +
      scale_shape_manual(values=shapes) +
      scale_color_manual(values=colors) +
      scale_y_continuous(statName, expand = expansion(mult = c(0.25,0.1))) +
      ggtitle(statName, subtitle=subtitle) +
      theme_replicates()
    return(p)
  }
  
  p_concBox <- plot_timesets("log_DNA_conc")
  p_tfBox <- plot_timesets("TF")
  p_fragBox <- plot_timesets("Fragment_size_ratio")
  
  # Assemble
  row1 <- plot_grid(pDist_conc, pDist_tf, pDist_frag, labels=LETTERS[1:3], nrow=1, label_size = 12)
  row2 <- plot_grid(p_concBox, p_tfBox, p_fragBox, labels=LETTERS[4:6], nrow=1, label_size = 12)
  
  final_grid <- plot_grid(row1, row2, labels = NULL, ncol=1, label_size = 12, rel_heights = c(1, 0.75))
  
  fname <- paste0(pdfdir, "Fig_S17_", filename_suffix, ".pdf")
  ggsave(final_grid, filename=fname, width=6.5, height=5)
  print(paste("Saved:", fname))
}

# --- Generate Figures ---
generate_stime_figure(res1$data_long, res1$data_box, "All")
generate_stime_figure(res2$data_long, res2$data_box, "Filtered")