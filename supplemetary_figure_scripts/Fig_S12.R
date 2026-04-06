library(tidyverse)
library(ggpubr)
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

# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

outputfile <- paste0(pdfdir, "Table.Fig_S12.stats.csv")
write_lines(c("Statistical results for analyses presented in Fig. S12"), outputfile, append=FALSE)

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

# --- Data Loading ---
raw <- as_tibble(read.csv(infile, header=T, sep="\t", na.strings = c("NA","N/A","","UNK"))) 
raw <- raw %>% filter(Sample_ID != "")
if("X" %in% colnames(raw)) raw$X <- NULL

# Recode Status
raw <- raw %>% mutate(
  Cancer_status = if_else(Cancer_type != "Healthy", "cancer", "healthy"),
  log_DNA_conc = log10(as.numeric(DNA_conc))
)

# Filter for Lymphoma & First Sample
pdAll <- raw %>% 
  filter(First_Sample) %>% 
  filter(Cohort != "ARM") %>%  # <--- EXCLUDE ARM COHORT
  filter(Cancer_type == "Lymphoma") %>%
  filter(!is.na(Stage)) %>%
  mutate(Stage = str_replace(Stage, "Stage ", "")) %>%
  # Keep only relevant stages
  filter(Stage %in% c("III", "IV", "V", "3", "4", "5")) %>%
  mutate(Stage = factor(Stage, levels = c("III", "IV", "V", "3", "4", "5"))) %>%
  # Normalize to Roman Numerals
  mutate(Stage = fct_recode(Stage, "III"="3", "IV"="4", "V"="5")) %>%
  # Drop unused levels
  mutate(Stage = droplevels(Stage)) %>%
  select(ID, Stage, log_DNA_conc, TF, Fragment_size_ratio) %>%
  pivot_longer(c(log_DNA_conc, TF, Fragment_size_ratio), names_to="stat", values_to="value") %>%
  filter(!is.na(value)) %>%
  left_join(rename_stats, by="stat")

# --- Define Comparisons ---
# Create explicit combinations
comps_list <- list(c("III", "IV"), c("IV", "V"), c("III", "V"))


# ==============================================================================
# STATISTICAL ANALYSIS (Pairwise: Wilcoxon, Permutation, Robust)
# ==============================================================================

calc_pairwise_stats <- function(data) {
  # Get all unique metrics
  metrics <- unique(data$stat)
  
  all_res <- map_dfr(metrics, function(m) {
    d_sub <- data %>% filter(stat == m)
    
    # 1. Base Wilcoxon test using rstatix
    # Explicitly use rstatix::wilcox_test to avoid conflict with coin::wilcox_test
    # Use p.adjust.method = "fdr" to generate p.adj column
    res <- d_sub %>% 
      rstatix::wilcox_test(value ~ Stage, p.adjust.method = "fdr") %>%
      mutate(stat = m)
    
    # 2. Add Permutation and Robust P-values
    # Iterate through the rows created by wilcox_test
    res$p_perm <- NA
    res$p_robust <- NA
    
    for(i in 1:nrow(res)) {
      g1 <- res$group1[i]
      g2 <- res$group2[i]
      
      # Subset data for this pair
      pair_data <- d_sub %>% 
        filter(Stage %in% c(g1, g2)) %>% 
        mutate(Stage = factor(Stage, levels=c(g1, g2)))
      
      # Check N
      if(nrow(pair_data) >= 4 && length(unique(pair_data$Stage)) == 2) {
        
        # Permutation (Independence Test)
        try({
          perm_obj <- coin::independence_test(value ~ Stage, data = pair_data, distribution="asymptotic")
          res$p_perm[i] <- as.numeric(coin::pvalue(perm_obj))
        }, silent=TRUE)
        
        # Robust (Yuen's test for independent samples)
        try({
          rob_obj <- WRS2::yuen(value ~ Stage, data = pair_data, tr=0.2)
          res$p_robust[i] <- rob_obj$p.value
        }, silent=TRUE)
      }
    }
    return(res)
  })
  
  return(all_res)
}

# Calculate Stats
stats_results <- calc_pairwise_stats(pdAll)

# Save to CSV
# Included p_fdr = p.adj
csv_out <- stats_results %>% 
  left_join(rename_stats %>% select(stat, statName), by="stat") %>%
  select(statName, group1, group2, n1, n2, p_wilcox=p, p_fdr=p.adj, p_perm, p_robust) %>% 
  mutate(across(starts_with("p_"), format_p_csv))

write.table(csv_out, file=outputfile, append=TRUE, sep=",", row.names=FALSE)


# ==============================================================================
# PLOTTING
# ==============================================================================

# Add y-positions for the brackets
stats_results <- stats_results %>% 
  left_join(pdAll %>% group_by(stat) %>% summarize(max_val = max(value, na.rm=T)), by="stat") %>%
  mutate(
    y.position = case_when(
      (group1 == "III" & group2 == "IV") ~ max_val * 1.05,
      (group1 == "IV" & group2 == "V")   ~ max_val * 1.15,
      (group1 == "III" & group2 == "V")  ~ max_val * 1.25,
      TRUE ~ max_val * 1.1
    )
  )

# Colors
stage_colors <- c("III" = "#ef3b2c", "IV" = "#a50f15", "V" = "#000000")

theme_stage <- function(){ 
  theme_cowplot(12) %+replace%
    theme(
      plot.title = element_text(hjust=0, size=8, face="bold", margin=margin(0,0,4,0)),
      plot.subtitle = element_text(hjust=0, size=7, margin=margin(0,0,4,0), lineheight=1),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=7, angle=90),
      axis.text.y = element_text(size=6, margin=margin(0,0,0,4)),
      axis.text.x = element_text(size=7, margin=margin(2,0,0,0))
    )
}

make_stage_plot <- function(metric) {
  pd <- pdAll %>% filter(stat == metric)
  stats_sub <- stats_results %>% filter(stat == metric)
  
  # Format p-values as numbers instead of stars
  stats_sub <- stats_sub %>% 
    mutate(p.adj.signif = case_when(
      p < 0.001 ~ formatC(p, format = "e", digits = 1),
      TRUE ~ formatC(p, format = "f", digits = 3)
    ))
  
  xlabels <- pd %>% 
    group_by(Stage) %>% 
    count() %>% 
    mutate(xlabel = paste0(Stage, "\n(N=", n, ")"))
  
  p <- ggplot(pd, aes(x=Stage, y=value)) +
    geom_boxplot(aes(fill=Stage), color="#525252", alpha=0.5, outlier.shape=NA, width=0.5) +
    geom_jitter(aes(color=Stage), alpha=0.75, width=0.1, height=0, size=1) +
    
    # Use Manual P-values
    stat_pvalue_manual(
      stats_sub, 
      label = "p.adj.signif", 
      tip.length = 0.01,
      size = 2.5 # Slightly smaller to fit text
    ) +
    
    ggtitle(unique(pd$statName), subtitle="Lymphoma by Stage") +
    scale_fill_manual(values = stage_colors) +
    scale_color_manual(values = stage_colors) +
    scale_x_discrete(labels = setNames(xlabels$xlabel, xlabels$Stage)) +
    theme_stage()
  
  # Y-axis scaling
  max_y <- max(pd$value, na.rm=T)
  
  if(metric == "log_DNA_conc") {
    p <- p + scale_y_continuous(unique(pd$axisLabel), breaks=c(-1:3), labels=10**c(-1:3), 
                                limits = c(min(pd$value), max_y * 1.3))
  } else if (metric == "TF") {
    p <- p + scale_y_continuous(unique(pd$axisLabel), breaks=c(0,0.5,1), 
                                limits = c(0, max(1.1, max_y * 1.3)))
  } else {
    p <- p + scale_y_continuous(unique(pd$axisLabel), 
                                limits = c(min(pd$value), max_y * 1.3))
  }
  
  return(p)
}

# Generate Plots
p_conc <- make_stage_plot("log_DNA_conc")
p_tf   <- make_stage_plot("TF")
p_frag <- make_stage_plot("Fragment_size_ratio")

# Assemble
pdfname <- paste0(pdfdir, "Fig_S12.pdf")
grid <- plot_grid(p_conc, p_tf, p_frag, labels="AUTO", ncol=3, label_size=12)
ggsave(grid, filename=pdfname, width=8, height=3.5)
print(paste("Made", pdfname))