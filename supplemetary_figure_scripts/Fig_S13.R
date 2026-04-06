library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(boot) # For CIs
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

outputfile <- paste0(pdfdir, "Table.Fig_S13.stats.csv")
write_lines(c("Statistical results for analyses presented in Fig. S13"), outputfile, append=FALSE)

rename_stats <- tibble(
  order=c(2,1,3),
  statName=c("Tumor fraction","cfDNA concentration","Fragment size ratio"),
  axisLabel=c("Tumor fraction","cfDNA concentration (ng/mL)","Fragment size ratio"),
  stat=c("TF","log_DNA_conc","Fragment_size_ratio") 
) %>% arrange(order)

# --- Helper: Formatter ---
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
raw <- as_tibble(read.csv(infile, header=T, sep="\t", na.strings = c("NA","N/A","","UNK"))) 
raw <- raw %>% filter(Sample_ID != "")
if("X" %in% colnames(raw)) raw$X <- NULL

# FILTER: Keep only NA, Rep1, or High replicates
if("replicate" %in% colnames(raw)) {
  raw <- raw %>% filter(is.na(replicate) | replicate %in% c("Rep1", "High"))
}

# Format Date
raw <- raw %>% mutate(Date_of_sample=str_replace(Date_of_sample," \\(approx\\)","")) %>% 
  mutate(Date_of_sample=parse_date_time(Date_of_sample,"m/d/y"))

# Recode Status & Log Transform
raw <- raw %>% mutate(
  Cancer_status = if_else(Cancer_type != "Healthy", "cancer", "healthy"),
  log_DNA_conc = log10(as.numeric(DNA_conc)),
  TF = as.numeric(TF),
  Fragment_size_ratio = as.numeric(Fragment_size_ratio),
  Tumor_size_cm = as.numeric(Tumor_size_cm)
)

# Propagate Tumor Size across ID
# Some samples might have tumor size recorded in one row but not the others
# We take the max tumor size recorded for the dog to ensure we capture the data for all samples
raw <- raw %>% 
  group_by(ID) %>% 
  mutate(Tumor_size_cm = ifelse(all(is.na(Tumor_size_cm)), NA, max(Tumor_size_cm, na.rm = TRUE))) %>%
  ungroup()

# --- Prepare Base Data for Analysis ---
# Cancer Only, Valid Tumor Size, Valid Metrics
base_data <- raw %>% 
  filter(Cancer_status == "cancer") %>% 
  filter(!is.na(Tumor_size_cm)) %>%
  # Create Grouping: Lymphoma vs Other
  mutate(Cancer_Group = if_else(Cancer_type == "Lymphoma", "Lymphoma", "Other cancers")) %>%
  select(ID, First_Sample, Cancer_Group, Tumor_size_cm, log_DNA_conc, TF, Fragment_size_ratio) %>%
  pivot_longer(c(log_DNA_conc, TF, Fragment_size_ratio), names_to="stat", values_to="value") %>%
  filter(!is.na(value)) %>%
  # Apply Log2 to Fragment Size Ratio for plotting consistency, protecting against <= 0
  mutate(plot_value = if_else(stat == "Fragment_size_ratio" & value > 0, log2(value), value)) %>%
  # Filter out -Inf or NaNs if any were produced (e.g. log of 0 or negative)
  filter(is.finite(plot_value)) %>%
  left_join(rename_stats, by="stat")

# Factor ordering for plot legend
base_data$Cancer_Group <- factor(base_data$Cancer_Group, levels=c("Lymphoma", "Other cancers"))

# --- Define Analysis Sets ---
pd_first <- base_data %>% filter(First_Sample)


# ==============================================================================
# ANALYSIS: FIRST SAMPLE ONLY (Spearman Correlation)
# ==============================================================================

# Boot function for Spearman
boot_spearman <- function(data, indices) {
  d <- data[indices,]
  return(cor(d$Tumor_size_cm, d$value, method = "spearman", use = "complete.obs"))
}

calc_group_stats_spearman <- function(data, group_name) {
  data %>%
    group_by(stat) %>%
    group_modify(~ {
      if(nrow(.x) < 3) return(tibble(n=nrow(.x), estimate=NA, p=NA, ci_low=NA, ci_high=NA))
      
      # Standard Test
      est <- cor.test(.x$Tumor_size_cm, .x$value, method="spearman", exact=FALSE)
      
      # Bootstrap CI
      set.seed(42)
      boot_res <- boot(data = .x, statistic = boot_spearman, R = 1000)
      ci <- tryCatch(boot.ci(boot_res, type="perc")$percent[4:5], error=function(e) c(NA,NA))
      
      tibble(
        Analysis = "First Sample (Spearman)",
        Group = group_name,
        n = nrow(.x),
        estimate = est$estimate,
        p = est$p.value,
        ci_low = ci[1],
        ci_high = ci[2]
      )
    })
}

# Run Spearman
stats_spearman_group <- pd_first %>% group_by(Cancer_Group) %>% group_modify(~ calc_group_stats_spearman(.x, as.character(.y$Cancer_Group)))
stats_spearman_all   <- calc_group_stats_spearman(pd_first, "All cancers")
final_stats <- bind_rows(stats_spearman_group, stats_spearman_all) %>%
  mutate(
    plot_label = paste0(Group, ": R=", round(estimate, 2), ", p=", format_p_plot(p))
  )


# ==============================================================================
# SAVE STATS & PLOTTING
# ==============================================================================

# Save Stats with Headers enabled
write.table(final_stats %>% 
              mutate(across(c(p, estimate, ci_low, ci_high), format_p_csv)), 
            file=outputfile, append=TRUE, sep=",", row.names=FALSE, col.names = TRUE)


# Define Colors
colors2 <- c("#ef3b2c", "#8073ac")

theme_scatter <- function(){ 
  theme_cowplot(12) %+replace%
    theme(
      plot.title = element_text(hjust=0, size=8, face="bold", margin=margin(0,0,4,0)),
      plot.subtitle = element_text(hjust=0, size=7, margin=margin(0,0,4,0), lineheight=1),
      legend.position = "none",
      axis.title = element_text(size=7),
      axis.text = element_text(size=6, margin=margin(2,0,0,4))
    )
}

# Plotting Function
make_scatter_plot <- function(data_source, stats_source, metric_name) {
  pd <- data_source %>% filter(stat == metric_name)
  stats_sub <- stats_source %>% filter(stat == metric_name & Group != "All cancers")
  
  # Y-limits
  max_y <- max(pd$plot_value, na.rm=T)
  min_y <- min(pd$plot_value, na.rm=T)
  range_y <- max_y - min_y
  
  # Label Placement
  stats_sub <- stats_sub %>% 
    mutate(
      y_pos = if_else(Group == "Lymphoma", max_y + (range_y * 0.15), max_y + (range_y * 0.05)),
      x_pos = min(pd$Tumor_size_cm, na.rm=T)
    )
  
  p <- ggplot(pd, aes(x=Tumor_size_cm, y=plot_value, group=Cancer_Group)) +
    geom_point(aes(color=Cancer_Group, shape=Cancer_Group), alpha=0.5, size=1.5) +
    # Use standard LM smoothing for visual trend
    geom_smooth(method="lm", aes(color=Cancer_Group, fill=Cancer_Group), alpha=0.15, size=0.5) +
    
    # Stats Labels
    geom_text(data=stats_sub, aes(x=x_pos, y=y_pos, label=plot_label, color=Group), 
              hjust=0, size=2.25, fontface="bold", show.legend=FALSE) +
    
    scale_color_manual(values=colors2) +
    scale_fill_manual(values=colors2) +
    scale_shape_manual(values=c(19, 17)) + 
    scale_x_continuous("Tumor size (cm)") +
    
    ggtitle(unique(pd$statName)) +
    theme_scatter()
  
  # Y-axis scaling (Matches metric types)
  if(metric_name == "log_DNA_conc") {
    p <- p + scale_y_continuous(unique(pd$axisLabel), breaks=c(-1:3), labels=10**c(-1:3), 
                                expand=expansion(mult=c(0.05, 0.25))) 
  } else if (metric_name == "Fragment_size_ratio") {
    breaks_frag <- c(0.25, 0.5, 1, 2, 4)
    p <- p + scale_y_continuous(unique(pd$axisLabel), 
                                breaks=log2(breaks_frag), labels=breaks_frag,
                                expand=expansion(mult=c(0.05, 0.25)))
  } else {
    p <- p + scale_y_continuous(unique(pd$axisLabel), 
                                expand=expansion(mult=c(0.05, 0.25)))
  }
  
  return(p)
}

# --- Row 1: First Sample (Spearman) ---
p1_conc <- make_scatter_plot(pd_first, final_stats, "log_DNA_conc")
p1_tf   <- make_scatter_plot(pd_first, final_stats, "TF")
p1_frag <- make_scatter_plot(pd_first, final_stats, "Fragment_size_ratio")

# Assemble
pdfname <- paste0(pdfdir, "Fig_S13.pdf")
grid <- plot_grid(p1_conc, p1_tf, p1_frag, ncol=3, labels="AUTO", label_size=12)

ggsave(grid, filename=pdfname, width=8, height=3) 
print(paste("Made", pdfname))