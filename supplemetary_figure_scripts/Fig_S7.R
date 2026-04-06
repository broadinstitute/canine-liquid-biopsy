library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(boot) # Added for CIs
# This package automatically finds the root folder of the downloaded project.
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ==============================================================================
# 1. SETUP & FILE PATHS
# ==============================================================================

# Define the input files
infile <- here("data", "BB_Supplementary_Data_1.txt")
depthfile <- here("data", "READDEPTH.txt")

# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

outputfile  <- paste0(pdfdir, "Table.Fig_S7.stats.csv")
write_lines(c("Statistical results for analyses presented in Fig. S7 (Spearman + Bootstrap CI)"), outputfile, append=FALSE)

# --- Helper: Scientific Notation Formatter (CSV) ---
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

# --- Helper: Plot Label Formatter ---
format_p_plot <- function(p) {
  sapply(p, function(val) {
    if (is.na(val)) return("NA")
    # Fix for 0.0e+00 - check for machine epsilon or literal zero
    if (val < 2.2e-16) return("< 2.2e-16")
    
    if (val < 0.001) {
      formatC(val, format = "e", digits = 1)
    } else {
      formatC(val, format = "f", digits = 3)
    }
  })
}

# --- Data Loading ---
# Columns: ID, cfDNA_TF, chr, start, end, tumor_log2_ratio, plasma_log2_ratio
ratios <- as_tibble(read.csv(depthfile, header=T, sep="\t")) %>% 
  rename(
    log2_tumor = tumor_log2_ratio,
    log2_cfDNA = plasma_log2_ratio
  ) %>%
  select(ID, log2_tumor, log2_cfDNA) %>% 
  filter(!is.na(log2_tumor) & !is.na(log2_cfDNA))


# ==============================================================================
# STATS: SPEARMAN + BOOTSTRAP CI (Per ID)
# ==============================================================================

# Bootstrap Function
boot_spearman <- function(data, indices) {
  d <- data[indices,]
  return(cor(d$log2_tumor, d$log2_cfDNA, method = "spearman", use = "complete.obs"))
}

stats_per_id <- ratios %>%
  group_by(ID) %>%
  group_modify(~ {
    # Check sample size (need at least 3 points for correlation)
    if(nrow(.x) < 3) return(tibble(n=nrow(.x), rho=NA, p=NA, ci_low=NA, ci_high=NA))
    
    # Standard Test
    est <- cor.test(.x$log2_tumor, .x$log2_cfDNA, method="spearman", exact=FALSE)
    
    # Bootstrap CI
    set.seed(42)
    boot_res <- boot(data = .x, statistic = boot_spearman, R = 1000)
    ci <- tryCatch(boot.ci(boot_res, type="perc")$percent[4:5], error=function(e) c(NA,NA))
    
    tibble(
      n = nrow(.x),
      rho = est$estimate,
      p = est$p.value,
      ci_low = ci[1],
      ci_high = ci[2]
    )
  }) %>%
  mutate(
    # String for Plot (R and P only)
    str_plot = paste0("R=", round(rho, 2), "; p=", format_p_plot(p)),
    # String for CSV (Full details)
    str_csv = paste0("Spearman: R=", round(rho, 2), " [", round(ci_low, 2), "-", round(ci_high, 2), "]")
  )

# Write detailed stats to CSV
csv_out <- stats_per_id %>%
  select(ID, n, rho, ci_low, ci_high, p, str_csv) %>%
  mutate(across(c(p, rho, ci_low, ci_high), format_p_csv))

write.table(csv_out, file=outputfile, append=TRUE, sep=",", row.names=FALSE)


# ==============================================================================
# PLOTTING
# ==============================================================================

# Prepare annotation dataframe (Merge stats with positioning coordinates)
pcorr <- ratios %>% 
  group_by(ID) %>% 
  summarize(
    log2_cfDNA = max(log2_cfDNA), 
    log2_tumor = min(log2_tumor)
  ) %>% 
  full_join(stats_per_id, by="ID")

### Set theme for these three plots
theme_compare <- function(){ 
  theme_cowplot(12) %+replace%
    theme(
      plot.title = element_text(hjust=0,size=8,face="bold",margin=margin(0,0,4,0)),
      plot.subtitle = element_text(hjust=0,size=7,margin=margin(0,0,4,0),lineheight = 1),
      axis.title = element_text(size=7),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = "grey70", linewidth = 0.25),
      axis.text = element_text(size=6,margin=margin(2,0,0,4)),
      legend.position=c(0.98,0.02),
      legend.justification = c("right", "bottom"),
      legend.margin = margin(3,3,3,3),
      legend.text = element_text(size=4.5,vjust=0.5),legend.key.size = unit(0.20, 'cm'),
      legend.title = element_blank(),
      legend.box.background = element_rect(color="#d9d9d9", size=0.5))
}

breaks <- c(-4,-2,-1,0,1,2,4)

p <- ggplot(ratios, aes(x=log2_tumor, y=log2_cfDNA)) 
p <- p + geom_abline(color="#525252", alpha=0.75, linetype=2, size=0.5)
p <- p + geom_point(color="#a50f15", shape=1, alpha=0.5)
p <- p + geom_text(aes(label=str_plot), data=pcorr, hjust=0, vjust=0, size=2, nudge_y=0.1)
p <- p + facet_wrap(~ID, scales="free", ncol=4)
p <- p + scale_x_continuous("copy number (tumor)", breaks=breaks, labels=2**breaks, expand = expansion(mult = c(0.15,0.15)))
p <- p + scale_y_continuous("copy number (cfDNA)", breaks=breaks, labels=2**breaks, expand = expansion(mult = c(0.15,0.25)))
p <- p + theme_compare()
p_ratios <- p 

### Make cowplot
pdfname <- paste0(pdfdir, "Fig_S7.pdf")
grid <- plot_grid(p_ratios, ncol=1, label_size = 12)
ggsave(grid, filename=pdfname, width=6.5, height=9)
print(paste("Made", pdfname))