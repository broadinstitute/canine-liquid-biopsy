library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(coin)
library(WRS2)
library(lme4)
library(lmerTest)
library(scales)
# This package automatically finds the root folder of the downloaded project.
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ==============================================================================
# 1. SETUP & FILE PATHS
# ==============================================================================

# Define the input files
infile <- here("data", "BB_Supplementary_Data_1.txt")
infile_ids <- here("data", "RETROSPECTIVE_METS.txt")
  
# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

outputfile <- paste0(pdfdir, "Table.Fig_S25_permuted.stats.csv")
write_lines(c("Statistical results for analyses presented in Fig. S25 (All Samples - Directional + Permuted LMM + FDR)"), outputfile, append=FALSE)

rename_stats <- tibble(
  order=c(2,1,3),
  statName=c("Tumor fraction","cfDNA concentration","Fragment size ratio"),
  axisLabel=c("Tumor fraction","cfDNA concentration (ng/mL)","Fragment size ratio"),
  stat=c("TF","log_DNA_conc","Fragment_size_ratio") 
) %>% arrange(order)

# --- Helper: Formatters ---
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
  sapply(p, function(val) {
    if (is.na(val)) return("NA")
    if (val < 0.001) {
      formatC(val, format = "e", digits = 1)
    } else {
      formatC(val, format = "f", digits = 3)
    }
  })
}

# --- Manual Van Elteren Test (Stratified Wilcoxon) ---
perm_van_elteren_manual <- function(value, group, block, direction="greater", B=10000) {
  
  calc_W_strat <- function(v, g, b) {
    blocks <- unique(b)
    total_W <- 0
    for(blk in blocks) {
      idx <- which(b == blk)
      vals <- v[idx]
      grps <- g[idx]
      if(length(unique(grps)) < 2) next
      x <- vals[grps ==levels(g)[2]] # Met
      y <- vals[grps ==levels(g)[1]] # Prior
      w <- wilcox.test(x, y, exact=FALSE)$statistic
      total_W <- total_W + w
    }
    return(total_W)
  }
  
  group <- factor(group) 
  W_obs <- calc_W_strat(value, group, block)
  
  set.seed(42)
  W_perms <- replicate(B, {
    # Shuffle group indices within blocks
    # Safer implementation using split/unsplit to maintain original order:
    perm_g <- unsplit(lapply(split(group, block), sample), block)
    calc_W_strat(value, perm_g, block)
  })
  
  if(direction == "greater") {
    p_val <- (sum(W_perms >= W_obs) + 1) / (B + 1)
  } else {
    p_val <- (sum(W_perms <= W_obs) + 1) / (B + 1)
  }
  return(p_val)
}

# --- Permuted Linear Mixed Model (1-Sided) ---
perm_lmm_1sided <- function(data, value_var, group_var, id_var, direction, B=2000) {
  # Use B=2000 for LMM as it is computationally heavier than Wilcoxon
  
  f <- as.formula(paste(value_var, "~", group_var, "+ (1|", id_var, ")"))
  
  # Fit Observed
  m_obs <- tryCatch(lmer(f, data=data), error=function(e) NULL)
  if(is.null(m_obs)) return(NA)
  
  # FIX: Access coefficient by index (Row 2 is Group effect) to avoid naming issues
  coefs_obs <- summary(m_obs)$coefficients
  if(nrow(coefs_obs) < 2) return(NA)
  t_obs <- coefs_obs[2, "t value"]
  
  set.seed(42)
  t_null <- numeric(B)
  
  # Prepare vectors for speed
  vec_grp <- data[[group_var]]
  vec_id  <- data[[id_var]]
  
  for(i in 1:B) {
    # Shuffle Group within ID
    perm_g <- unsplit(lapply(split(vec_grp, vec_id), sample), vec_id)
    d_perm <- data
    d_perm[[group_var]] <- perm_g
    
    # Fit Permuted
    m_perm <- suppressMessages(suppressWarnings(tryCatch(lmer(f, data=d_perm), error=function(e) NULL)))
    
    if(!is.null(m_perm)) {
      coefs <- summary(m_perm)$coefficients
      if(nrow(coefs) >= 2) {
        t_null[i] <- coefs[2, "t value"]
      } else {
        t_null[i] <- NA
      }
    } else {
      t_null[i] <- NA
    }
  }
  
  # Filter NAs (failed convergence)
  t_null <- t_null[!is.na(t_null)]
  if(length(t_null) < 100) return(NA) # Return NA if too many perms failed
  
  # Calculate 1-sided P
  # NOTE: Direction logic flipped
  if(direction == "greater") {
    p_val <- (sum(t_null <= t_obs) + 1) / (length(t_null) + 1)
  } else {
    p_val <- (sum(t_null >= t_obs) + 1) / (length(t_null) + 1)
  }
  return(p_val)
}


# --- Data Loading ---
raw <- as_tibble(read.csv(infile, header=T, sep="\t", na.strings = c("NA","N/A","","UNK"))) 
raw <- raw %>% filter(Sample_ID != "")

raw <- raw %>% mutate(
  log_DNA_conc = log10(as.numeric(DNA_conc)),
  TF = as.numeric(TF),
  Fragment_size_ratio = as.numeric(Fragment_size_ratio)
)

# Load ID Mapping
meta_ids <- as_tibble(read.delim(infile_ids, header=T, sep="\t", stringsAsFactors=FALSE))

# FIX: Correct ID
meta_ids <- meta_ids %>% 
  mutate(across(everything(), ~str_replace(., "T-OS-0061-cfDNA--WK8", "T-OS-0061-cfDNA-1ng-WK8"))) %>%
  mutate(across(everything(), trimws))

# Reshape Mapping
prior_map <- meta_ids %>% 
  select(ID, starts_with("Previous")) %>% 
  pivot_longer(cols = starts_with("Previous"), values_to = "Sample_ID") %>% 
  filter(!is.na(Sample_ID) & Sample_ID != "") %>% 
  select(ID, Sample_ID) %>% 
  mutate(Group = "No evidence of disease")

met_map <- meta_ids %>% 
  select(ID, Sample_ID) %>% 
  filter(!is.na(Sample_ID) & Sample_ID != "") %>% 
  mutate(Group = "Retrospectively identified metastasis")

map_final <- bind_rows(prior_map, met_map)

# Merge
pdAll <- raw %>% 
  inner_join(map_final, by = "Sample_ID") %>% 
  mutate(ID = coalesce(ID.x, ID.y)) %>% 
  select(ID, Sample_ID, Group, replicate, ULP_input_DNA, log_DNA_conc, TF, Fragment_size_ratio)

# --- DATA PROCESSING ---
# Keep ALL samples (no averaging)
pdStats <- pdAll %>%
  pivot_longer(c(log_DNA_conc, TF, Fragment_size_ratio), names_to="stat", values_to="value") %>%
  filter(!is.na(value)) %>%
  left_join(rename_stats, by="stat")

pdStats$Group <- factor(pdStats$Group, levels = c("No evidence of disease", "Retrospectively identified metastasis"))


# ==============================================================================
# STATISTICAL ANALYSIS
# ==============================================================================

calc_retro_stats <- function(data) {
  
  res <- data %>%
    group_by(statName) %>%
    group_modify(~ {
      d_sub <- .x
      
      # Direction
      direction <- if_else(unique(d_sub$stat) == "Fragment_size_ratio", "less", "greater")
      
      # 1. Standard Linear Mixed Model (1-sided)
      m <- tryCatch(lmer(value ~ Group + (1|ID), data = d_sub), error=function(e) NULL)
      p_lmm_1sided <- NA
      
      if(!is.null(m)) {
        coefs <- summary(m)$coefficients
        
        # FIX: Access Row 2 (Variable effect) directly to avoid naming mismatch
        if(nrow(coefs) >= 2) {
          t_stat <- coefs[2, "t value"]
          
          # Check if df is available (lmerTest)
          if("df" %in% colnames(coefs)) {
            df <- coefs[2, "df"]
            # NOTE: Direction logic flipped
            if (direction == "greater") {
              p_lmm_1sided <- pt(t_stat, df, lower.tail = TRUE)
            } else {
              p_lmm_1sided <- pt(t_stat, df, lower.tail = FALSE)
            }
          } else {
            # Fallback to Normal approximation if df estimation fails
            # (Warning: less accurate for small N, but prevents NA)
            if (direction == "greater") {
              p_lmm_1sided <- pnorm(t_stat, lower.tail = TRUE)
            } else {
              p_lmm_1sided <- pnorm(t_stat, lower.tail = FALSE)
            }
          }
        }
      }
      
      # 2. Permuted Linear Mixed Model (1-sided)
      p_lmm_perm <- perm_lmm_1sided(d_sub, "value", "Group", "ID", direction, B=2000)
      
      # 3. Manual Van Elteren Test (Permuted, 1-sided)
      p_van_elteren <- perm_van_elteren_manual(d_sub$value, d_sub$Group, d_sub$ID, direction, B=10000)
      
      tibble(
        N_obs = nrow(d_sub),
        N_dogs = n_distinct(d_sub$ID),
        direction = direction,
        p_lmm_1sided = p_lmm_1sided,
        p_lmm_perm = p_lmm_perm,
        p_van_elteren = p_van_elteren
      )
    }) %>%
    ungroup() %>%
    # Add FDR Adjustment across the 3 metrics
    mutate(
      p_lmm_fdr = p.adjust(p_lmm_1sided, method = "fdr"),
      p_lmm_perm_fdr = p.adjust(p_lmm_perm, method = "fdr"),
      p_van_elteren_fdr = p.adjust(p_van_elteren, method = "fdr")
    )
  
  return(res)
}

print("Running Statistical Analysis (LMM + Permutations)...")
stats_results <- calc_retro_stats(pdStats)

# Save to CSV
write_lines(c("", "Comparison: No evidence vs Retrospective Met (All Samples - Directional Tests with FDR)"), outputfile, append=TRUE)
write.table(stats_results %>% mutate(across(starts_with("p_"), format_p_csv)), 
            file=outputfile, append=TRUE, sep=",", row.names=FALSE, col.names=TRUE)


# ==============================================================================
# PLOTTING
# ==============================================================================

# Join stats for labels
pdPlot <- pdStats %>% left_join(stats_results, by="statName") %>%
  mutate(
    subtitle = paste0("LMM (1-sided) p = ", format_p_plot(p_lmm_1sided),
                      "\nLMM (Permuted) p = ", format_p_plot(p_lmm_perm),
                      "\nVan Elteren p = ", format_p_plot(p_van_elteren))
  )

# Calculate Means for Connecting Lines
pdMeans <- pdPlot %>%
  group_by(ID, stat, Group) %>%
  summarize(mean_value = mean(value, na.rm=TRUE), .groups="drop")

colors_retro <- c("No evidence of disease" = "#4d4d4d", "Retrospectively identified metastasis" = "#a50f15")

theme_retro <- function(){ 
  theme_cowplot(12) %+replace%
    theme(
      plot.title = element_text(hjust=0, size=8, face="bold", margin=margin(0,0,4,0)),
      plot.subtitle = element_text(hjust=0, size=6, margin=margin(0,0,4,0), lineheight=1.2),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=7, angle=90),
      axis.text.y = element_text(size=6, margin=margin(0,0,0,4)),
      axis.text.x = element_text(size=7, margin=margin(2,0,0,0), angle=0, hjust=0.5)
    )
}

make_retro_plot <- function(metric) {
  pd_dat <- pdPlot %>% filter(stat == metric)
  pd_lines <- pdMeans %>% filter(stat == metric)
  sub_txt <- unique(pd_dat$subtitle)
  
  p <- ggplot(pd_dat, aes(x=Group, y=value)) +
    # Lines (connect means)
    geom_line(data=pd_lines, aes(x=Group, y=mean_value, group=ID), color="#bdbdbd", alpha=0.6, size=0.5) +
    
    # Boxplots
    geom_boxplot(aes(fill=Group), alpha=0.5, outlier.shape=NA, width=0.5, color="#525252") +
    
    # Points (All samples)
    geom_jitter(aes(color=Group), width=0.1, height=0, size=1.5, alpha=0.8) +
    
    scale_fill_manual(values = colors_retro) +
    scale_color_manual(values = colors_retro) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) + 
    ggtitle(unique(pd_dat$statName), subtitle = sub_txt) +
    theme_retro()
  
  if(metric == "log_DNA_conc") {
    p <- p + scale_y_continuous(unique(pd_dat$axisLabel), breaks=c(-1:3), labels=10**c(-1:3), expand=expansion(mult=c(0.1, 0.1)))
  } else if (metric == "TF") {
    p <- p + scale_y_continuous(unique(pd_dat$axisLabel), expand=expansion(mult=c(0.1, 0.15)), breaks=scales::pretty_breaks())
  } else {
    p <- p + scale_y_continuous(unique(pd_dat$axisLabel), expand=expansion(mult=c(0.1, 0.1)))
  }
  
  return(p)
}

p_conc <- make_retro_plot("log_DNA_conc")
p_tf   <- make_retro_plot("TF")
p_frag <- make_retro_plot("Fragment_size_ratio")

pdfname <- paste0(pdfdir, "Fig_S25.pdf")
grid <- plot_grid(p_tf, p_frag, labels="AUTO", ncol=2, label_size=12)
ggsave(grid, filename=pdfname, width=6.5, height=4) # Slightly taller for 3 lines of subtitle
print(paste("Made", pdfname))