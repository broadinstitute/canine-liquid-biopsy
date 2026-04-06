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
infile <- here("data", "CONMUTS.txt")
mathfile <- here("data", "CONMATH.txt")
namefile <- here("data", "BB_Supplementary_Data_1.txt")

# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

outputfile <- paste0(pdfdir, "Table.Fig_S9.stats.csv")
write_lines(c("Statistical results for analyses presented in Fig. S9"), outputfile, append=FALSE)

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

# --- Helper: Plot Label Formatter ---
format_p_plot <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ formatC(p, format = "e", digits = 2),
    TRUE ~ formatC(p, format = "f", digits = 3)
  )
}

# --- Data Loading ---

# 1. Dog Names
dognames <- as_tibble(read.csv(namefile, header=T, sep="\t", na.strings = c("NA","N/A","","UNK","?"))) %>% 
  select(ID) %>% 
  distinct() %>% 
  mutate(dogname = ID) 

# 2. Math Data
math <- as_tibble(read.csv(mathfile, header=T, sep="\t")) 

# 3. Mutation Data
raw <- as_tibble(read.csv(infile, header=T, sep="\t", na.strings=c("","NA","NI","UNK"), strip.white=TRUE, stringsAsFactors=TRUE)) %>% 
  filter(FILTER=="PASS") %>% 
  select(-FILTER)

# Rename annotation columns gracefully
if("ANN.0..EFFECT" %in% colnames(raw)) {
  raw <- raw %>% rename(EFFECT=ANN.0..EFFECT, IMPACT=ANN.0..IMPACT, GENE=ANN.0..GENE, GENEID=ANN.0..GENEID, BIOTYPE=ANN.0..BIOTYPE, HGVS_C=ANN.0..HGVS_C, HGVS_P=ANN.0..HGVS_P)
}

# Identify unique mutations
mutations <- raw %>% select(CHROM, POS, REF, ALT, GENE, GENEID) %>% distinct()
locations <- mutations %>% 
  arrange(CHROM, POS) %>% 
  mutate(mutID=row_number())

mutations <- mutations %>% left_join(locations) 
raw <- raw %>% inner_join(mutations)

# --- METRIC CALCULATION ---

# 1. Count-based metrics (Total, Unique, Others)
counts_long <- raw %>%
  group_by(ID, Sample_type) %>%
  summarize(total_count = n(), .groups="drop")

# Calculate overlap (nboth)
muts_distinct <- raw %>% select(ID, Sample_type, CHROM, POS, REF, ALT) %>% distinct()
overlap_counts <- muts_distinct %>%
  group_by(ID, CHROM, POS, REF, ALT) %>%
  summarize(n_types = n_distinct(Sample_type), .groups="drop") %>%
  filter(n_types > 1) %>%
  group_by(ID) %>%
  summarize(nboth = n(), .groups="drop")

# Merge counts and overlap
counts_wide <- counts_long %>%
  pivot_wider(names_from = Sample_type, values_from = total_count, values_fill = 0) %>%
  left_join(overlap_counts, by="ID") %>%
  mutate(nboth = replace_na(nboth, 0))

# Calculate Metrics
metrics_counts <- counts_wide %>%
  pivot_longer(c(tumor, cfDNA), names_to = "type", values_to = "total_count") %>%
  mutate(
    total.detected = total_count,
    unique.detected = total_count - nboth,
    others.detected = if_else(total_count > 0, nboth / total_count, 0)
  ) %>%
  select(ID, type, total.detected, unique.detected, others.detected) %>%
  pivot_longer(cols = c(total.detected, unique.detected, others.detected), names_to = "name", values_to = "value")

# 2. Value-based metrics (MATH, Median.AF)
if("Sample_type" %in% colnames(math)) {
  metrics_math <- math %>%
    select(ID, Sample_type, any_of(c("MATH", "Median.AF"))) %>%
    rename(type = Sample_type) %>%
    pivot_longer(any_of(c("MATH", "Median.AF")), names_to = "name", values_to = "value")
} else {
  metrics_math <- tibble(ID=character(), type=character(), name=character(), value=numeric())
}

# Combine all metrics
pdAll <- bind_rows(metrics_counts, metrics_math) %>%
  filter(!is.na(value)) %>%
  filter(type %in% c("tumor", "cfDNA"))

pdAll$type <- factor(pdAll$type, levels=c("tumor", "cfDNA"))


# ==============================================================================
# STATISTICAL ANALYSIS (T-test, Permutation, Robust)
# ==============================================================================

paired_stats <- pdAll %>%
  group_by(name) %>%
  group_modify(~ {
    d_wide <- .x %>% 
      select(ID, type, value) %>%
      pivot_wider(names_from=type, values_from=value) %>%
      drop_na(tumor, cfDNA)
    
    if(nrow(d_wide) < 3) {
      return(tibble(n=nrow(d_wide), p_ttest=NA, p_perm=NA, p_robust=NA))
    }
    
    # 1. Paired T-test (Used for Plot)
    p_ttest <- tryCatch(t.test(d_wide$tumor, d_wide$cfDNA, paired=TRUE)$p.value, error=function(e) NA)
    
    # 2. Permutation (Validation)
    d_long <- d_wide %>% pivot_longer(c(tumor, cfDNA), names_to="type", values_to="val") %>%
      mutate(type = factor(type, levels=c("tumor", "cfDNA")), ID = factor(ID))
    
    p_perm <- tryCatch({
      as.numeric(pvalue(symmetry_test(val ~ type | ID, data=d_long, distribution="asymptotic")))
    }, error=function(e) NA)
    
    # 3. Robust (Validation)
    p_robust <- tryCatch({
      WRS2::yuend(d_wide$tumor, d_wide$cfDNA, tr=0.2)$p.value
    }, error=function(e) NA)
    
    tibble(
      n = nrow(d_wide),
      p_ttest = p_ttest,
      p_perm = p_perm,
      p_robust = p_robust
    )
  }) %>%
  ungroup() %>%
  # FDR Adjustment for T-test P-values
  mutate(p_ttest_adj = p.adjust(p_ttest, method = "fdr"))

# Write Stats (p_ttest_adj will be included via starts_with("p_"))
write.table(paired_stats %>% mutate(across(starts_with("p_"), format_p_csv)), 
            file=outputfile, append=TRUE, sep=",", row.names=FALSE)


# ==============================================================================
# PLOTTING
# ==============================================================================

# Prepare Data
pdPlot <- pdAll %>% 
  left_join(dognames, by="ID") %>%
  mutate(dogname = ifelse(is.na(dogname), ID, dogname))

# Connecting Lines
linesPlot <- pdPlot %>%
  select(ID, dogname, name, type, value) %>%
  pivot_wider(names_from=type, values_from=value) %>%
  drop_na(tumor, cfDNA) %>%
  mutate(type="tumor", type2="cfDNA", value=tumor, value2=cfDNA) 

# Titles Mapping
titles <- tibble(
  name = c("Median.AF","MATH","total.detected","unique.detected","others.detected"),
  ytitle = c("median allele frequency","MATH score","# mutations","# mutations","fraction mutations"),
  plottitle = c("Heterogeneity:\nMedian allele frequency","Heterogeneity:\nMATH score","Number of unique\nmutations detected","Total number of\nmutations detected","Fraction of other\nmutations detected")
) %>% 
  left_join(paired_stats, by="name") %>%
  mutate(
    # Using T-test P-value for plot
    subtitle = paste("Paired t-test p = ", format_p_plot(p_ttest), "\n", sep="")
  )

theme_paired <- function(){ 
  theme_cowplot(12) %+replace%
    theme(
      plot.title = element_text(hjust=0,size=7,face="bold"),
      plot.subtitle = element_text(size=7,hjust=0),
      axis.title.x=element_blank(),
      axis.title.y = element_text(vjust=1,hjust=0.5,size=6,angle=90),
      axis.text = element_text(size=6),
      legend.position = "NONE")
}

make_paired_plot <- function(label) {
  dat <- pdPlot %>% filter(name == label)
  lin <- linesPlot %>% filter(name == label)
  tit <- titles %>% filter(name == label)
  
  if(nrow(dat) == 0) return(ggplot() + theme_void())
  
  p <- ggplot(dat, aes(x=type, y=value)) + 
    geom_boxplot(outlier.shape = NA, width=0.5, alpha=0.25, color="#525252", fill="#969696", linewidth=0.5) +
    geom_text(aes(label=dogname, color=dogname), x=2.35, data=dat %>% filter(type=="cfDNA"), hjust=0, size=1.75) +
    geom_segment(aes(x=type, xend=type2, y=value, yend=value2, color=dogname), data=lin, linewidth=0.5) +
    geom_point(aes(fill=dogname), color="#525252", shape=21) +
    ggtitle(tit$plottitle, subtitle=tit$subtitle) +
    scale_y_continuous(tit$ytitle) +
    scale_x_discrete("", expand = expansion(mult = c(0.75, 1.75))) +
    theme_paired()
  
  return(p)
}

# Generate Plots
p_total  <- make_paired_plot("total.detected")
p_unique <- make_paired_plot("unique.detected")
p_over   <- make_paired_plot("others.detected")
p_AF     <- make_paired_plot("Median.AF")
p_math   <- make_paired_plot("MATH")

# Assemble
pdfname <- paste0(pdfdir, "Fig_S9.pdf")
row2 <- plot_grid(p_total, p_unique, p_over, p_AF, p_math, ncol=3, label_size = 12, labels=LETTERS)
grid <- plot_grid(row2, ncol=1, label_size = 12)

ggsave(grid, filename=pdfname, width=6, height=6)
print(paste("Made", pdfname))

# Save Descriptive Stats
desc_stats <- pdPlot %>% 
  group_by(type, name) %>% 
  summarize(n=n(), min=min(value, na.rm=T), max=max(value, na.rm=T), median=median(value, na.rm=T), mean=mean(value, na.rm=T))

write_lines(c("", "Descriptive Statistics"), outputfile, append=TRUE)
write.table(desc_stats, file=outputfile, append=TRUE, sep=",", row.names=FALSE)