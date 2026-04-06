library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(flextable)
library(coin) 
library(WRS2)
library(boot) # Added for CIs
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

# Setup
rename_stats <- tibble(
  order=c(2,1,3),
  statName=c("Tumor fraction","cfDNA concentration","Fragment size ratio"),
  axisLabel=c("Tumor fraction","cfDNA concentration (ng/mL)","Fragment size ratio"),
  stat=c("TF","log_DNA_conc","Fragment_size_ratio") 
) %>% arrange(order)

outputfile <- paste0(pdfdir, "Fig_S4_stats.csv")
write_lines(c("Statistical results for analyses presented in Fig. S4"),outputfile,append=FALSE)

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

# --- Data Loading ---
raw <- as_tibble(read.csv(infile,header=T,sep='\t',na.strings = c("NA","N/A","","UNK"))) 

# Create Cancer_status from Cancer_type
raw <- raw %>% mutate(
  Cancer_status = if_else(Cancer_type != "Healthy", "cancer", "healthy"),
  Cancer_status2 = Cancer_status
)

# Ensure numeric types
raw <- raw %>% mutate(
  log_DNA_conc = log10(as.numeric(DNA_conc)),
  TF = as.numeric(TF),
  Fragment_size_ratio = as.numeric(Fragment_size_ratio)
) %>% select(-DNA_conc)


# ==============================================================================
# DATA PREPARATION
# ==============================================================================

# 1. Select relevant columns (Include Plasma_ID to separate aliquots)
pdAll <- raw %>% 
  filter(!is.na(replicate)) %>%
  select(ID, Plasma_ID, Cancer_status, replicate, log_DNA_conc, TF, Fragment_size_ratio) %>%
  distinct() %>% 
  pivot_longer(c(log_DNA_conc, TF, Fragment_size_ratio), names_to = "stat", values_to = "value") %>%
  filter(!is.na(value)) 

# 2. Pivot to Wide
pdAll <- pdAll %>% 
  pivot_wider(
    names_from = replicate, 
    values_from = value,
    values_fn = function(x) x[1]
  )

# 3. Filter for valid pairs
if("Rep1" %in% colnames(pdAll) & "Rep2" %in% colnames(pdAll)) {
  pdAll <- pdAll %>% filter(!is.na(Rep1) & !is.na(Rep2))
} else {
  stop("Columns Rep1 and Rep2 not created. Check input data for 'replicate' column.")
}

# 4. Add metadata back
pdAll <- pdAll %>% mutate(shapeID=if_else(Cancer_status=="healthy","healthy dogs",ID)) 
pdAll <- pdAll %>% left_join(rename_stats, by = "stat")


# ==============================================================================
# HELPER: BOOTSTRAP FUNCTION
# ==============================================================================
boot_spearman <- function(data, indices) {
  d <- data[indices,]
  return(cor(d$Rep1, d$Rep2, method = "spearman", use = "complete.obs"))
}

# ==============================================================================
# PART 1: OVERALL SPEARMAN (Label: "Spearman")
# ==============================================================================

stats_overall <- pdAll %>% 
  ungroup() %>% 
  group_by(statName) %>% 
  group_modify(~ {
    # Standard Test
    est <- cor.test(.x$Rep1, .x$Rep2, method="spearman", exact=FALSE)
    
    # Bootstrap CI
    set.seed(42)
    boot_res <- boot(data = .x, statistic = boot_spearman, R = 1000)
    boot_ci <- boot.ci(boot_res, type = "perc")
    
    tibble(
      n = nrow(.x),
      p = est$p.value,
      R = est$estimate,
      ci_lower = boot_ci$percent[4],
      ci_upper = boot_ci$percent[5]
    )
  }) %>%
  mutate(
    # Detailed string for CSV (Optional, but good for checking)
    str_csv = paste0("Spearman: R=", round(R, 2), " [", round(ci_lower, 2), "-", round(ci_upper, 2), "]"),
    # Simple string for Plot
    str_overall = paste("Spearman: R=",round(R,2),", p=",format(p,digits=2),sep="")
  )

write_lines(c("",paste(LETTERS[1],"Overall Spearman correlation analysis (All Replicates)")),outputfile,append=TRUE)
write.table(stats_overall %>% select(-str_overall), file=outputfile, append=TRUE, sep=",", row.names=FALSE)


# ==============================================================================
# PART 2: CANCER CASES SPEARMAN (Label: "Partial Spearman (cancer)")
# ==============================================================================

stats_cancer <- pdAll %>% 
  filter(Cancer_status == "cancer") %>%
  ungroup() %>% 
  group_by(statName) %>% 
  group_modify(~ {
    # Standard Test
    est <- cor.test(.x$Rep1, .x$Rep2, method="spearman", exact=FALSE)
    
    # Bootstrap CI
    set.seed(42)
    boot_res <- boot(data = .x, statistic = boot_spearman, R = 1000)
    boot_ci <- boot.ci(boot_res, type = "perc")
    
    tibble(
      n_cases = nrow(.x),
      p_cases = est$p.value,
      R_cases = est$estimate,
      ci_lower_cases = boot_ci$percent[4],
      ci_upper_cases = boot_ci$percent[5]
    )
  }) %>%
  mutate(
    # Simple string for Plot
    str_cancer = paste("Partial Spearman (cancer): R=",round(R_cases,2),", p=",format(p_cases,digits=2),sep="")
  )

write_lines(c("",paste(LETTERS[2],"Spearman correlation analysis (Cancer Cases Only)")),outputfile,append=TRUE)
write.table(stats_cancer %>% select(-str_cancer), file=outputfile, append=TRUE, sep=",", row.names=FALSE)

# Join stats for plotting labels
final_stats <- stats_overall %>% left_join(stats_cancer, by="statName")


# ==============================================================================
# PLOTTING
# ==============================================================================

theme_replicates <- function(){ 
  theme_cowplot(12) %+replace%
    theme(
      plot.title = element_text(hjust=0,size=8,face="bold",margin=margin(0,0,4,0)),
      plot.subtitle = element_text(hjust=0,size=6,margin=margin(0,0,4,0),lineheight = 1.1),
      plot.margin = unit(c(0.2, 0.3, 0.1, 0.1), "cm"),
      axis.title = element_text(size=7),
      axis.text = element_text(size=6,margin=margin(2,0,0,4)),
      legend.position=c(0.98,0.02),
      legend.justification = c("right", "bottom"),
      legend.margin = margin(3,3,3,3),
      legend.text = element_text(size=4.5,vjust=0.5),legend.key.size = unit(0.20, 'cm'),
      legend.title = element_blank(),
      legend.box.background = element_rect(color="#d9d9d9", size=0.5))
}

# Helper to generate plots
make_rep_plot <- function(data, metric, stats_df) {
  pd <- data %>% filter(stat == metric)
  stat_name <- unique(pd$statName)
  
  # Get labels
  row <- stats_df %>% filter(statName == stat_name)
  subtitle <- paste(row$str_overall, row$str_cancer, sep="\n")
  
  # Dynamic shape handling
  n_shapes <- length(unique(pd$shapeID))
  shape_vals <- if(n_shapes > 5) rep(16, n_shapes) else c(0,1,2,5,16)
  
  p <- ggplot(pd, aes(x=Rep1, y=Rep2)) +
    geom_smooth(method = "lm", color="#525252", size=0.5, linewidth=0.25, alpha=0.25) +
    geom_point(aes(color=shapeID, shape=shapeID), alpha=0.75, size=1.25) +
    scale_color_manual(values=c("#cb181d","#cb181d","#cb181d","#cb181d","#252525")) +
    scale_fill_manual(values=c("#cb181d","#cb181d","#cb181d","#cb181d","#252525")) +
    scale_shape_manual(values=shape_vals) +
    scale_x_continuous("Replicate 1") +
    scale_y_continuous("Replicate 2") +
    ggtitle(stat_name, subtitle=subtitle) +
    theme_replicates()
  
  return(p)
}

# Generate Plots
p_compare_conc <- make_rep_plot(pdAll, "log_DNA_conc", final_stats)
p_compare_tf   <- make_rep_plot(pdAll, "TF", final_stats)
p_compare_frag <- make_rep_plot(pdAll, "Fragment_size_ratio", final_stats)

# Assemble
row_compare <- plot_grid(p_compare_conc, p_compare_tf, p_compare_frag, labels=LETTERS[1:3], ncol=3, label_size=12)
pdfname <- paste0(pdfdir, "Fig_S4.pdf")
grid <- plot_grid(row_compare, labels=NULL, ncol=1, label_size=12) 

ggsave(grid, filename=pdfname, width=6.5, height=2.5)
print(paste("Made", pdfname))