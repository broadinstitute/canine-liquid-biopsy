# Load libraries
library(survival)
library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)
library(scales) # For label formatting
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

# ==========================================
# 1. Load and Preprocess Data
# ==========================================

# Read data
data <- read.delim(infile, stringsAsFactors = FALSE, check.names = FALSE)

# Fix duplicate column names
names(data) <- make.unique(names(data))

# -------------------------------------------------------
# COLUMN NAME CLEANUP
# -------------------------------------------------------
tumor_col_idx <- grep("^tumor_size_cm$", names(data), ignore.case = TRUE)
if (length(tumor_col_idx) > 0) {
  names(data)[tumor_col_idx[1]] <- "tumor_size_cm"
}

# Clean whitespace
data$Cancer_type <- trimws(data$Cancer_type)
data$Sample_time_point <- trimws(data$Sample_time_point)

# RECODING: Map negative time points and windowing
data$Sample_time_point[grepl("Day\\s*-[0-9]", data$Sample_time_point)] <- "Day 0"
data$Sample_time_point[data$Sample_time_point %in% c("Day 8", "Day 10")] <- "Day 7"

# Filter: Cancer_type == "Lymphoma"
df_analysis_base <- data %>%
  filter(Cancer_type == "Lymphoma")

print(paste("Total Lymphoma Dogs Found:", nrow(df_analysis_base)))

# Survival: 1 = Event
df_analysis_base$PFS_censor <- as.numeric(df_analysis_base$PFS_censor)
df_analysis_base$Event <- df_analysis_base$PFS_censor 
df_analysis_base$PFS <- as.numeric(df_analysis_base$PFS)

# -------------------------------------------------------
# FORMAT COVARIATES & TRANSFORMATIONS
# -------------------------------------------------------
df_analysis_base <- df_analysis_base %>%
  mutate(
    Age = as.numeric(Age_at_sample),
    Weight = as.numeric(Weight_kg),
    Sex = as.factor(Sex),
    Disease_Status = as.factor(Disease_status),
    Stage = as.factor(Stage),  
    BCS = as.numeric(BCS),
    Tumor_Size = as.numeric(tumor_size_cm),
    Breed = as.factor(Breed),
    Log10_DNA_conc = log10(as.numeric(DNA_conc))
  )

# ==========================================
# 2. Variable Definitions & Helpers
# ==========================================

metric_cols <- c("Log10_DNA_conc", "TF", "Fragment_size_ratio")
metric_cols <- intersect(metric_cols, names(df_analysis_base))
covariates <- c("Age", "Sex", "Weight", "Disease_Status", "Stage", "BCS", "Tumor_Size","Breed")

# Helper: Permutation Test (Univariate Robustness)
# Uses LIKELIHOOD RATIO TEST (LRT) STATISTIC
run_permutation_uni_clean <- function(time, event, val, n_perm = 1000, metric_name="Metric") {
  
  # 1. Create a local dataframe
  local_data <- data.frame(Time = time, Event = event, Value = val)
  
  # 2. Fit Observed Model & Extract LRT Statistic (2 * diff(loglik))
  obs_stat <- tryCatch({
    fit <- coxph(Surv(Time, Event) ~ Value, data = local_data)
    if (length(fit$loglik) < 2) {
      NA # Return NA to obs_stat, don't exit function
    } else {
      2 * (fit$loglik[2] - fit$loglik[1]) # Return value to obs_stat
    }
  }, error = function(e) {
    message(paste("Permutation Obs Error:", e$message))
    NA
  })
  
  if (is.na(obs_stat)) return(NA)
  
  # 3. Run Permutations
  perm_stats <- numeric(n_perm)
  valid_perms <- 0
  
  for(i in 1:n_perm) {
    local_data$Shuffled <- sample(local_data$Value)
    
    # Use suppressWarnings to avoid separation warnings spam
    try_res <- try({
      fit_perm <- suppressWarnings(coxph(Surv(Time, Event) ~ Shuffled, data = local_data))
      if (length(fit_perm$loglik) >= 2) {
        2 * (fit_perm$loglik[2] - fit_perm$loglik[1])
      } else {
        NA
      }
    }, silent=TRUE)
    
    if (!inherits(try_res, "try-error") && !is.na(try_res)) {
      perm_stats[i] <- try_res
      valid_perms <- valid_perms + 1
    } else {
      perm_stats[i] <- NA
    }
  }
  
  # 4. Calculate Empirical P-value
  if (valid_perms < 10) {
    message("Permutation Failed: Too few valid permutations.")
    return(NA)
  }
  
  valid_stats <- perm_stats[!is.na(perm_stats)]
  
  # Count how many random permutations had a stronger LRT than observed
  n_exceed <- sum(valid_stats >= obs_stat)
  
  # Pseudo-count correction (Standard P-value calculation)
  pval <- (n_exceed + 1) / (length(valid_stats) + 1)
  
  # Sanity cap
  if (pval > 1) pval <- 1
  
  return(pval)
}

# Helper: Permutation Test (Multivariable Robustness)
# Uses Wald Statistic (Z-score)
run_permutation_mv_clean <- function(f, data_in, var_name = "Value", n_perm = 1000) {
  
  # 1. Fit Observed Model
  obs_fit <- try(coxph(f, data = data_in), silent=TRUE)
  if (inherits(obs_fit, "try-error")) return(NA)
  
  coef_summ <- summary(obs_fit)$coefficients
  if (!(var_name %in% rownames(coef_summ))) return(NA)
  obs_stat <- abs(coef_summ[var_name, "z"])
  
  # 2. Run Permutations
  perm_stats <- numeric(n_perm)
  valid_perms <- 0
  perm_data <- data_in
  
  for(i in 1:n_perm) {
    perm_data[[var_name]] <- sample(data_in[[var_name]])
    
    try_res <- try({
      fit_perm <- suppressWarnings(coxph(f, data = perm_data))
      summ_perm <- summary(fit_perm)$coefficients
      if (var_name %in% rownames(summ_perm)) {
        abs(summ_perm[var_name, "z"])
      } else {
        NA
      }
    }, silent=TRUE)
    
    if (!inherits(try_res, "try-error") && !is.na(try_res) && try_res < 20) {
      perm_stats[i] <- try_res
      valid_perms <- valid_perms + 1
    } else {
      perm_stats[i] <- NA
    }
  }
  
  # 3. Calculate P-value
  if (valid_perms < 10) return(NA)
  valid_stats <- perm_stats[!is.na(perm_stats)]
  
  n_exceed <- sum(valid_stats >= obs_stat)
  pval <- (n_exceed + 1) / (length(valid_stats) + 1)
  
  if (pval > 1) pval <- 1
  return(pval)
}

# Helper: Bootstrap Test (Multivariable Robustness)
run_bootstrap_mv <- function(f, data_in, n_boot = 500) {
  res_hr <- numeric(n_boot)
  n <- nrow(data_in)
  
  for(i in 1:n_boot) {
    indices <- sample(1:n, n, replace=TRUE)
    boot_dat <- data_in[indices, ]
    
    try({
      fit <- suppressWarnings(coxph(f, data=boot_dat))
      if(any(abs(coef(fit)) > 15, na.rm=TRUE)) {
        res_hr[i] <- NA
      } else {
        coefs <- coef(fit)
        if("Value" %in% names(coefs)) {
          res_hr[i] <- exp(coefs[["Value"]])
        } else {
          res_hr[i] <- NA
        }
      }
    }, silent=TRUE)
  }
  
  res_hr <- res_hr[!is.na(res_hr)]
  if(length(res_hr) < 50) return(c(NA, NA)) 
  
  return(quantile(res_hr, probs=c(0.025, 0.975)))
}

# ==========================================
# 3. ANALYSIS LOOP
# ==========================================

results_df <- data.frame()
time_points <- c("Day 0", "Day 7") 

print("Starting Analysis...")

for (tp in time_points) {
  
  print(paste("--- PROCESSING:", tp, "---"))
  df_time <- df_analysis_base %>% filter(Sample_time_point == tp)
  
  if(nrow(df_time) == 0) { next }
  
  # A. Univariate Models (With Permutation Check)
  for (metric in metric_cols) {
    sub_uni <- na.omit(df_time %>% select(PFS, Event, all_of(metric)) %>% rename(Value = all_of(metric)))
    if(nrow(sub_uni) > 5) {
      tryCatch({
        # Standard Cox Model
        fit_uni <- coxph(Surv(PFS, Event) ~ Value, data = sub_uni)
        stats_uni <- tidy(fit_uni, exponentiate = TRUE, conf.int = TRUE)
        
        # Permutation Test (Using LRT statistic)
        perm_p <- run_permutation_uni_clean(sub_uni$PFS, sub_uni$Event, sub_uni$Value, n_perm=1000, metric_name=metric)
        
        results_df <- rbind(results_df, data.frame(
          Time_Point = tp,
          Model_Type = "Univariate",
          Metric = metric,
          HR = stats_uni$estimate,
          Lower_CI = stats_uni$conf.low,
          Upper_CI = stats_uni$conf.high,
          P_value = stats_uni$p.value,
          Permutation_P = perm_p,
          Boot_Lower_CI = NA, # N/A for Univariate
          Boot_Upper_CI = NA,
          N_Dogs = fit_uni$n,
          N_Events = fit_uni$nevent
        ))
      }, error = function(e) { message(paste("Uni Error", metric, ":", e$message))})
    }
  }
  
  # B. Multivariable Models (Fully Adjusted + Bootstrap + Permutation)
  for (metric in metric_cols) {
    sub_multi <- na.omit(df_time %>% select(PFS, Event, all_of(metric), all_of(covariates)) %>% rename(Value = all_of(metric)))
    if(nrow(sub_multi) > 10) {
      tryCatch({
        valid_covs <- c()
        for (cov in covariates) {
          if (is.numeric(sub_multi[[cov]]) && var(sub_multi[[cov]], na.rm=TRUE) == 0) next
          if (length(unique(sub_multi[[cov]])) < 2) next
          if ((is.factor(sub_multi[[cov]]) || is.character(sub_multi[[cov]])) && any(tapply(sub_multi$Event, sub_multi[[cov]], sum) == 0, na.rm=TRUE)) next
          valid_covs <- c(valid_covs, cov)
        }
        
        if (length(valid_covs) > 0) {
          f <- as.formula(paste("Surv(PFS, Event) ~ Value +", paste(valid_covs, collapse=" + ")))
          fit_multi <- coxph(f, data = sub_multi)
          stats_multi <- tidy(fit_multi, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "Value")
          
          if(nrow(stats_multi) > 0) {
            # Run Bootstrap
            boot_cis <- run_bootstrap_mv(f, sub_multi, n_boot=500)
            
            # Run Permutation (MV)
            perm_p_mv <- run_permutation_mv_clean(f, sub_multi, var_name="Value", n_perm=1000)
            
            results_df <- rbind(results_df, data.frame(
              Time_Point = tp,
              Model_Type = "Multivariable",
              Metric = metric,
              HR = stats_multi$estimate,
              Lower_CI = stats_multi$conf.low,
              Upper_CI = stats_multi$conf.high,
              P_value = stats_multi$p.value,
              Permutation_P = perm_p_mv,
              Boot_Lower_CI = boot_cis[1],
              Boot_Upper_CI = boot_cis[2],
              N_Dogs = fit_multi$n,
              N_Events = fit_multi$nevent
            ))
          }
        }
      }, error = function(e) { message(paste("Multi Error", metric, ":", e$message)) })
    }
  }
}

# ==========================================
# 4. Save Statistics
# ==========================================
if (nrow(results_df) > 0) {
  out_file <- paste0(pdfdir, "Fig_S26_stats.txt")
  
  # Format the dataframe for saving
  # Prevent scientific notation in P-values
  save_df <- results_df
  save_df$P_value <- sprintf("%.5f", save_df$P_value)
  save_df$Permutation_P <- sprintf("%.5f", save_df$Permutation_P)
  
  write.table(save_df, out_file, sep="\t", quote=FALSE, row.names=FALSE)
  print(paste("Statistics saved to:", out_file))
}

# ==========================================
# 5. Generate Plots (Combined)
# ==========================================

if (nrow(results_df) > 0) {
  
  results_df$P_fmt <- sprintf("%.3g", results_df$P_value)
  
  results_df$Label <- paste0("N=", results_df$N_Dogs, 
                             ", Ev=", results_df$N_Events, 
                             ", p=", results_df$P_fmt)
  
  results_df$Metric <- factor(results_df$Metric, levels = rev(metric_cols))
  results_df$Model_Type <- factor(results_df$Model_Type, levels = c("Univariate", "Multivariable"))
  results_df$Time_Point <- factor(results_df$Time_Point, levels = c("Day 0", "Day 7"))
  
  p <- ggplot(results_df, aes(x = HR, y = Metric, color = Model_Type, group = Model_Type)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.3, position = position_dodge(width = 0.7)) +
    geom_point(size = 3, position = position_dodge(width = 0.7)) +
    geom_text(aes(label = Label, x = Upper_CI), hjust = -0.1, vjust = 0.5, size = 3, position = position_dodge(width = 0.7), show.legend = FALSE) +
    scale_color_manual(values = c("Univariate" = "#1f77b4", "Multivariable" = "#d62728")) +
    scale_x_log10(labels = label_number(drop0trailing = TRUE)) + 
    facet_wrap(~Time_Point, ncol = 2, scales = "free_x") +
    coord_cartesian(clip = "off") +
    scale_x_continuous(trans='log10', expand = expansion(mult = c(0.05, 0.4))) +
    
    theme_bw() +
    labs(
      title = "Cox PH Models: Day 0 vs Day 7",
      subtitle = "Comparison of Biomarker Metrics (Univariate & Multivariable)",
      x = "Hazard Ratio (Log Scale)", 
      y = "",
      color = "Model"
    ) +
    theme(
      legend.position = "bottom",
      plot.margin = unit(c(1, 1, 1, 1), "lines"),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 12, face = "bold")
    )
  
  print(p)
  
  filename_png <- paste0(pdfdir,"S26.png")
  ggsave(filename_png, plot = p, width = 14, height = 7, dpi = 300)
  print(paste("Saved Combined Plot (PNG):", filename_png))
  
  filename_pdf <- paste0(pdfdir,"S26.pdf")
  ggsave(filename_pdf, plot = p, width = 14, height = 7)
  print(paste("Saved Combined Plot (PDF):", filename_pdf))
}