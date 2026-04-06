library(tidyverse)
library(cowplot)
library(ggpubr)
library(rstatix)
library(lubridate)
library(lme4)      # For Mixed Models
library(lmerTest)  # For p-values in Mixed Models
library(coin)      # For Permutation Tests
library(robustbase) # For Robust Regression
# This package automatically finds the root folder of the downloaded project.
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ==============================================================================
# 1. SETUP & FILE PATHS
# ==============================================================================

# Define the input files
infile <- here("data", "BB_Supplementary_Data_1.txt")
infile_cbc <- here("data", "SCBC.txt")

# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

pdfname  <- paste(pdfdir,"Fig_7_LMM_AllSamples.pdf",sep="")
outputfile <- paste0(pdfdir, "Table.Fig_7_LMM_AllSamples.stats.csv")

write_lines(c("Statistical results for analyses presented in Fig. 7 (Longitudinal LMM - All Samples)"), outputfile, append=FALSE)

rename_stats <- tibble(
  order    = c(2,1,3),
  statistic= c("Tumor fraction","cfDNA (ng/mL)","Fragment size ratio"),
  stat     = c("TF","log_DNA_conc","Fragment_size_ratio")
) %>%
  arrange(order) %>%
  select(-order)

cohorts <- c("OS_Longitudinal",
             "LSA_U01_treatment_cohort_2",
             "LSA_U01_treatment_cohort_1")

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
  # Force numeric to prevent "unsupported type" error in formatC
  p <- as.numeric(p) 
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ formatC(p, format = "e", digits = 1),
    TRUE ~ formatC(p, format = "f", digits = 3)
  )
}

# --- Data Loading ---
raw <- as_tibble(read.csv(infile, header = TRUE, sep = "\t", na.strings = c("NA","N/A","","UNK")))

if (!"Short_ID" %in% colnames(raw)) {
  raw <- raw %>% mutate(Short_ID = NA_character_)
}

raw <- raw %>% filter(Cohort %in% cohorts)

raw <- raw %>%
  mutate(
    Date_of_sample  = str_replace(Date_of_sample," \\(approx\\)",""),
    Date_of_sample  = parse_date_time(Date_of_sample,"m/d/y")
  )

raw <- raw %>%
  mutate(
    Date_of_diagnosis = str_remove(Date_of_diagnosis,"10/28/2011 and "),
    Date_of_diagnosis = str_replace(Date_of_diagnosis," \\(approx\\)",""),
    Date_of_diagnosis = parse_date_time(Date_of_diagnosis,"m/d/y")
  )

raw <- raw %>%
  mutate(
    Date_of_progressive_disease2    = parse_date_time(Date_of_progressive_disease,"m/d/y"),
    Last_known_date_if_still_alive = parse_date_time(Last_known_date_if_still_alive,"m/d/y")
  )

raw <- raw %>% mutate(log_DNA_conc = log10(DNA_conc))
raw <- raw %>% mutate(Cancer_status = if_else(Cancer_type!="Healthy","cancer","healthy"))

# convert Disease status in readable form and align terminology
# UPDATED: Ensured "no evidence of disease" is used
Disease_status <- tibble(
  Disease_status      = c("Before tx","disease present","CR","PR","SD","NED","Relapse","PD"),
  Disease_status_new = c("Before tx","disease present","no evidence of disease",
                         "partial resp.","stable disease","no evidence of disease",
                         "relapse","progressive disease")
)

Disease_status <- raw %>%
  select(Disease_status) %>%
  distinct() %>%
  inner_join(Disease_status)

Disease_status$status_order <- c(4,5,3,1)

raw <- raw %>%
  left_join(Disease_status) %>%
  rename(Disease_status_short = Disease_status) %>%
  rename(Disease_status       = Disease_status_new)

Disease_status <- raw %>%
  select(Disease_status,status_order) %>%
  distinct()

## make everything referenced on Date of first sample
raw <- raw %>%
  group_by(ID) %>%
  summarize(Date_of_first_sample = min(Date_of_sample)) %>%
  full_join(raw)

raw <- raw %>%
  mutate(timepoint = as.Date(Date_of_sample) - as.Date(Date_of_first_sample))

# Filter duplicates to keep max DNA input (same as before)
inputdna <- raw %>%
  filter(!is.na(ULP_input_DNA)) %>%
  select(ID,Sample_ID,Plasma_ID,ULP_input_DNA) %>%
  group_by(ID,Plasma_ID) %>%
  summarize(max = max(ULP_input_DNA)) %>%
  inner_join(raw, relationship = "many-to-many")

inputdna <- inputdna %>%
  ungroup() %>%
  mutate(keep = if_else(ULP_input_DNA==max,TRUE,FALSE)) %>%
  select(ID,Sample_ID,keep)

raw <- raw %>%
  left_join(inputdna, relationship = "many-to-many") %>%
  filter(is.na(ULP_input_DNA) | keep) %>%
  select(-keep)

dogs <- raw %>%
  select(ID,Cancer_type,Short_ID) %>%
  distinct() %>%
  mutate(longname = if_else(is.na(Short_ID) | Short_ID==ID, ID, paste(Short_ID,"\n",ID,sep="")))

raw <- raw %>%
  select(ID,Sample_ID,log_DNA_conc,TF,Fragment_size_ratio,timepoint,Disease_status) %>%
  distinct()

# Identify Relapse
progressive <- raw %>%
  filter(Disease_status=="progressive disease") %>%
  select(ID,timepoint) %>%
  group_by(ID) %>%
  summarize(timepoint = min(timepoint)) %>%
  mutate(relapse = TRUE)

raw <- raw %>%
  left_join(progressive) %>%
  replace_na(list(relapse = FALSE))

#### PREPARE DATA ####
pdAll <- raw %>%
  select(ID,Sample_ID,timepoint,log_DNA_conc,TF,Fragment_size_ratio,Disease_status,relapse) %>%
  pivot_longer(c(log_DNA_conc,TF,Fragment_size_ratio)) %>%
  rename(stat = name)

pdAll <- dogs %>% select(ID,Cancer_type,longname) %>% inner_join(pdAll)

# Logic: Timepoint 0 is "disease present", otherwise use status
pdAll <- pdAll %>%
  mutate(plotbin = if_else(timepoint==0,"disease present",Disease_status))

pdAll <- pdAll %>% left_join(rename_stats)

## remove partial response
pdAll <- pdAll %>% filter(Disease_status!="partial resp.")

# UPDATED: Labels set to "no evidence of disease" without newlines for rotation
plotbins    <- c("disease present","no evidence of disease","progressive disease")
pdAll$plotbin <- factor(pdAll$plotbin, levels = plotbins)
xlabels <- c("disease present","no evidence of disease","progressive disease")

# ==============================================================================
# STATISTICAL ANALYSIS (Linear Mixed Models + Robust + Permutation)
# ==============================================================================
# Since we have multiple samples per dog per bin, we use LMM
# Model: Value ~ Status + (1|ID)

calc_lmm_stats <- function(data) {
  
  pairs_to_test <- list(
    c("disease present", "no evidence of disease"),
    c("no evidence of disease", "progressive disease")
  )
  
  res <- data %>%
    group_by(Cancer_type, stat) %>%
    group_modify(~ {
      
      curr_res <- map_dfr(pairs_to_test, function(pair) {
        g1 <- pair[1]
        g2 <- pair[2]
        
        # Filter for IDs that have BOTH states
        ids_g1 <- .x %>% filter(plotbin == g1) %>% pull(ID) %>% unique()
        ids_g2 <- .x %>% filter(plotbin == g2) %>% pull(ID) %>% unique()
        common_ids <- intersect(ids_g1, ids_g2)
        
        # FILTER: Ensure values are finite (No -Inf for logs)
        d_sub <- .x %>% 
          filter(plotbin %in% c(g1, g2)) %>%
          filter(ID %in% common_ids) %>%
          filter(is.finite(value)) %>% # CRITICAL FIX
          mutate(plotbin = factor(plotbin, levels = c(g1, g2))) 
        
        if(length(common_ids) < 3 || nrow(d_sub) < 4) {
          return(tibble(Group1=g1, Group2=g2, N_Dogs=length(common_ids), N_Samples=nrow(d_sub), 
                        p_lmm=NA, p_perm=NA, p_robust=NA))
        }
        
        # 1. Linear Mixed Model (Primary)
        # Fixed effect: Disease Status, Random effect: Dog ID
        p_lmm <- tryCatch({
          m <- lmer(value ~ plotbin + (1|ID), data = d_sub)
          coefs <- summary(m)$coefficients
          coefs[2, "Pr(>|t|)"]
        }, error = function(e) NA)
        
        # 2. Permutation Test (coin)
        # Explicitly use coin::approximate and handle factor levels
        p_perm <- tryCatch({
          d_perm <- d_sub %>% 
            mutate(
              ID = factor(ID), 
              plotbin = droplevels(factor(plotbin)) # Ensure clean factors
            ) %>% 
            as.data.frame()
          
          # independence_test is robust to unbalanced data/distributions
          it <- coin::independence_test(value ~ plotbin | ID, data = d_perm, 
                                        distribution = coin::approximate(nresample = 5000))
          as.numeric(coin::pvalue(it))
        }, error = function(e) NA)
        
        # 3. Robust Linear Regression (Robust Check)
        p_robust <- tryCatch({
          m_rob <- robustbase::lmrob(value ~ plotbin, data = d_sub, setting = "KS2014")
          summary(m_rob)$coefficients[2, "Pr(>|t|)"]
        }, error = function(e) NA)
        
        tibble(
          Group1 = g1, 
          Group2 = g2, 
          N_Dogs = length(common_ids), 
          N_Samples = nrow(d_sub),
          p_lmm = as.numeric(p_lmm),
          p_perm = as.numeric(p_perm),
          p_robust = as.numeric(p_robust)
        )
      })
      return(curr_res)
    }) %>%
    ungroup() %>%
    # Add FDR Adjustment (on the primary LMM p-value)
    mutate(p_adj = p.adjust(p_lmm, method = "fdr"))
  
  return(res)
}

# Run Stats
stats_out <- calc_lmm_stats(pdAll)

# Save to CSV (Suppress append warning)
write_lines(c("", "Longitudinal Comparison (LMM, Permutation, Robust, FDR)"), outputfile, append=TRUE)
suppressWarnings(
  write.table(stats_out %>% mutate(across(starts_with("p_"), format_p_csv)), 
              file=outputfile, append=TRUE, sep=",", row.names=FALSE)
)

# ==============================================================================
# PLOTTING
# ==============================================================================

theme_boxes <- function(){ 
  theme_cowplot(12) %+replace%
    theme(
      plot.title    = element_text(hjust=0,size=7,face="bold",
                                   margin=margin(0,0,4,0)),
      plot.subtitle = element_text(hjust=0,size=6,
                                   margin=margin(0,0,4,0),lineheight = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=8,angle=90,vjust=1),
      # UPDATED: Rotated 45 degrees, slightly larger size
      axis.text.x  = element_text(size=6, angle=45, hjust=1, vjust=1,
                                  margin=margin(2,0,0,4)),
      axis.text.y  = element_text(size=5,
                                  margin=margin(2,0,0,4)),
      strip.text    = element_text(size=6,hjust=0.5,
                                   face="bold",color="#f7f7f7"),
      strip.background = element_rect(fill="#525252"),
      panel.border     = element_rect(color = "#525252",
                                      fill  = NA,
                                      size  = 0.5),
      legend.position  = "none"
    )
}

## Helper: Get Dog Ns for X-axis labels
get_ns <- function(statIn) {
  pd <- pdAll %>% filter(stat == statIn, !is.na(value))
  ns <- pd %>%
    group_by(plotbin) %>%
    summarise(N = n_distinct(ID), .groups = "drop")
  
  lab_df <- tibble(
    plotbin    = factor(plotbins, levels = plotbins),
    base_label = xlabels
  ) %>%
    left_join(ns, by = "plotbin") %>%
    mutate(N = if_else(is.na(N), 0L, N))
  
  lab_df %>%
    mutate(label = paste0(base_label, "\nN=", N)) %>%
    pull(label)
}

# --- PLOT FUNCTION ---
plot_lmm <- function(statIn) {
  pd <- pdAll %>% filter(stat == statIn) %>% filter(!is.na(value))
  
  # Prepare Line Data: Averages per Dog per Bin
  # (Lines connect centroids to keep visual clutter down, points show variance)
  pdLines <- pd %>%
    group_by(ID, Cancer_type, plotbin) %>%
    summarize(value = mean(value, na.rm=TRUE), .groups="drop")
  
  # Identify IDs valid for specific comparisons for linking lines
  compare_set1 <- pdLines %>%
    filter(plotbin %in% plotbins[1:2]) %>%
    group_by(ID) %>%
    count() %>%
    filter(n==2) %>%
    select(ID)
  
  compare_set2 <- pdLines %>%
    filter(plotbin %in% plotbins[2:3]) %>%
    group_by(ID) %>%
    count() %>%
    filter(n==2) %>%
    select(ID)
  
  okIds <- compare_set1 %>% bind_rows(compare_set2) %>% distinct()
  pdLines_filtered <- pdLines %>% inner_join(okIds)
  
  # Prepare P-values for plotting
  # We extract them from stats_out and create a dataframe for geom_text
  pvals <- stats_out %>% 
    filter(stat == statIn) %>%
    # Filter NAs first to avoid formatting errors
    filter(!is.na(p_lmm)) %>%
    mutate(
      label = paste0("p=", format_p_plot(p_lmm)),
      x_pos = if_else(Group1 == "disease present", 1.5, 2.5),
      # Adjust Y position based on data max
      y_pos = Inf
    )
  
  # Get Max Y for placing p-values (approximate)
  max_y_vals <- pd %>% group_by(Cancer_type) %>% summarize(max_val = max(value, na.rm=T))
  pvals <- pvals %>% left_join(max_y_vals, by="Cancer_type")
  
  # Plotting
  p <- ggplot(pd, aes(x=plotbin, y=value))
  
  # 1. Raw Points (Jittered slightly to show multiple samples)
  p <- p + geom_jitter(aes(color=Cancer_type), shape=16, size=1, alpha=0.5, width=0.1, height=0)
  
  # 2. Connecting Lines (Using Averaged Data to show Trajectory)
  p <- p + geom_line(data=pdLines_filtered, aes(group=ID), color="#707070", alpha=0.4, size=0.4)
  
  # 3. P-values (Manual Text)
  # UPDATED: Smaller text size (1.8) and adjusted vjust (-1.0) for spacing
  p <- p + geom_text(data=pvals, aes(x=x_pos, y=max_val, label=label), 
                     vjust=-1.0, size=1.8, inherit.aes=FALSE)
  
  p <- p + scale_color_manual(values=c("#ef3b2c","#8073ac"))
  p <- p + facet_wrap(~Cancer_type, nrow=1, scales="free")
  
  labels_with_n <- get_ns(statIn)
  p <- p + scale_x_discrete("", breaks = xbreaks, labels = labels_with_n)
  
  # UPDATED: Decreased expansion to 0.15 (from 0.25)
  if (statIn=="log_DNA_conc") {
    p <- p + scale_y_continuous(unique(pd$statistic), breaks = c(0:4), labels = 10**c(0:4), expand = expansion(mult = c(0.05, 0.15)))
  } else {
    p <- p + scale_y_continuous(unique(pd$statistic), expand = expansion(mult = c(0.05, 0.15)))
  }
  
  p <- p + ggtitle(unique(pd$statistic), subtitle = "Longitudinal LMM (All Samples)")
  p <- p + theme_boxes()
  p
}

p_conc2 <- plot_lmm("log_DNA_conc")
p_TF2   <- plot_lmm("TF")
p_frag2 <- plot_lmm("Fragment_size_ratio")

row1 <- plot_grid(p_conc2, p_TF2, nrow=1, labels=LETTERS[1:2], rel_widths = c(1,2))
row2 <- plot_grid(p_frag2, nrow=1, ncol=2, labels=LETTERS[3], rel_widths = c(2,1))
grid <- plot_grid(row1, row2, ncol=1)

ggsave(grid, filename = pdfname, width = 4, height = 5)
print(paste("Made", pdfname))

# ==============================================================================
# BOXPLOT VERSION (SEPARATE PDF)
# ==============================================================================

plot_lmm_boxplot <- function(statIn) {
  pd <- pdAll %>% filter(stat == statIn) %>% filter(!is.na(value))
  
  # Prepare Line Data: Averages per Dog per Bin (same as above)
  pdLines <- pd %>%
    group_by(ID, Cancer_type, plotbin) %>%
    summarize(value = mean(value, na.rm=TRUE), .groups="drop")
  
  compare_set1 <- pdLines %>%
    filter(plotbin %in% plotbins[1:2]) %>%
    group_by(ID) %>%
    count() %>%
    filter(n==2) %>%
    select(ID)
  
  compare_set2 <- pdLines %>%
    filter(plotbin %in% plotbins[2:3]) %>%
    group_by(ID) %>%
    count() %>%
    filter(n==2) %>%
    select(ID)
  
  okIds <- compare_set1 %>% bind_rows(compare_set2) %>% distinct()
  pdLines_filtered <- pdLines %>% inner_join(okIds)
  
  # Prepare P-values
  pvals <- stats_out %>% 
    filter(stat == statIn) %>%
    # Filter NAs first
    filter(!is.na(p_lmm)) %>%
    mutate(
      label = paste0("p=", format_p_plot(p_lmm)),
      x_pos = if_else(Group1 == "disease present", 1.5, 2.5),
      y_pos = Inf
    )
  
  max_y_vals <- pd %>% group_by(Cancer_type) %>% summarize(max_val = max(value, na.rm=T))
  pvals <- pvals %>% left_join(max_y_vals, by="Cancer_type")
  
  # Plotting
  p <- ggplot(pd, aes(x=plotbin, y=value))
  
  # 1. Boxplot (Added)
  p <- p + geom_boxplot(aes(fill=Cancer_type), alpha=0.3, outlier.shape=NA, width=0.6)
  
  # 2. Connecting Lines (Fainter for boxplot version)
  p <- p + geom_line(data=pdLines_filtered, aes(group=ID), color="#707070", alpha=0.2, size=0.3)
  
  # 3. Jittered Points (Added)
  p <- p + geom_jitter(aes(color=Cancer_type), shape=16, size=0.8, alpha=0.6, width=0.15, height=0)
  
  # 4. P-values
  p <- p + geom_text(data=pvals, aes(x=x_pos, y=max_val, label=label), 
                     vjust=-1.0, size=1.8, inherit.aes=FALSE)
  
  p <- p + scale_color_manual(values=c("#ef3b2c","#8073ac"))
  p <- p + scale_fill_manual(values=c("#ef3b2c","#8073ac"))
  p <- p + facet_wrap(~Cancer_type, nrow=1, scales="free")
  
  labels_with_n <- get_ns(statIn)
  p <- p + scale_x_discrete("", breaks = xbreaks, labels = labels_with_n)
  
  if (statIn=="log_DNA_conc") {
    p <- p + scale_y_continuous(unique(pd$statistic), breaks = c(0:4), labels = 10**c(0:4), expand = expansion(mult = c(0.05, 0.15)))
  } else {
    p <- p + scale_y_continuous(unique(pd$statistic), expand = expansion(mult = c(0.05, 0.15)))
  }
  
  p <- p + ggtitle(unique(pd$statistic), subtitle = "Longitudinal LMM (Boxplots)")
  p <- p + theme_boxes()
  p
}

p_conc_box <- plot_lmm_boxplot("log_DNA_conc")
p_TF_box   <- plot_lmm_boxplot("TF")
p_frag_box <- plot_lmm_boxplot("Fragment_size_ratio")

row1_box <- plot_grid(p_conc_box, p_TF_box, nrow=1, labels=LETTERS[1:2], rel_widths = c(1,2))
row2_box <- plot_grid(p_frag_box, nrow=1, ncol=2, labels=LETTERS[3], rel_widths = c(2,1))
grid_box <- plot_grid(row1_box, row2_box, ncol=1)

pdfname_box <- paste(pdfdir,"Fig_7_LMM_AllSamples_Boxplot.pdf",sep="")
ggsave(grid_box, filename = pdfname_box, width = 4, height = 5)
print(paste("Made", pdfname_box))