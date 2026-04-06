library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(flextable)
library(coin) 
library(WRS2)
library(boot) 
if (!requireNamespace("robustbase", quietly = TRUE)) install.packages("robustbase")
library(robustbase) 
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

outputfile <- paste0(pdfdir, "Fig_S5_stats.csv")
write_lines(c("Statistical results for analyses presented in Fig. S5"), outputfile, append=FALSE)

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
raw <- as_tibble(read.csv(infile, header=T, sep='\t', na.strings = c("NA","N/A","","UNK"))) 
raw <- raw %>% mutate(
  Cancer_status = if_else(Cancer_type != "Healthy", "cancer", "healthy"),
  Cancer_status2 = Cancer_status,
  log_DNA_conc = log10(as.numeric(DNA_conc)),
  TF = as.numeric(TF),
  Fragment_size_ratio = as.numeric(Fragment_size_ratio)
) %>% select(-DNA_conc)

# Filter for rows with ULP Input DNA data
raw_long <- raw %>% 
  filter(!is.na(ULP_input_DNA)) %>% 
  select(ID, Plasma_ID, Sample_ID, ULP_input_DNA, TF, Fragment_size_ratio) %>%
  group_by(Plasma_ID) %>% 
  filter(n() > 1) %>%
  ungroup() %>%
  pivot_longer(c(TF, Fragment_size_ratio), names_to = "stat", values_to = "value")

# --- Create Pairs (For Plotting & Spearman) ---
pairs <- raw_long %>% 
  select(-Sample_ID) %>% 
  distinct() 

pairs <- pairs %>% 
  rename(ULP_input_DNA_2 = ULP_input_DNA, value_2 = value) %>% 
  inner_join(pairs, by = c("ID", "Plasma_ID", "stat"), relationship = "many-to-many") %>%
  filter(ULP_input_DNA_2 > ULP_input_DNA) %>% 
  mutate(diff = ULP_input_DNA_2 - ULP_input_DNA) %>%
  rename(lower = value, higher = value_2)


# ==============================================================================
# STATS PART 1: SPEARMAN CORRELATION + BOOTSTRAP CI
# ==============================================================================

# Function for boot library
boot_spearman <- function(data, indices) {
  d <- data[indices,]
  return(cor(d$lower, d$higher, method = "spearman", use = "complete.obs"))
}

cor_stats <- pairs %>% 
  group_by(stat) %>% 
  group_modify(~ {
    # 1. Standard Spearman
    est <- cor.test(.x$lower, .x$higher, method = "spearman", exact=FALSE)
    
    # 2. Bootstrap 95% CI (Robustness check)
    set.seed(42)
    boot_res <- boot(data = .x, statistic = boot_spearman, R = 1000)
    boot_ci <- boot.ci(boot_res, type = "perc")
    
    tibble(
      n = nrow(.x),
      spearman_rho = est$estimate,
      p_spearman = est$p.value,
      ci_lower = boot_ci$percent[4],
      ci_upper = boot_ci$percent[5]
    )
  }) %>%
  # Create two label strings: one for CSV (detailed), one for Plot (simple)
  mutate(
    str_cor_csv = paste0("Spearman: R=", round(spearman_rho, 2), 
                         " [", round(ci_lower, 2), "-", round(ci_upper, 2), "]",
                         "; p=", format(p_spearman, digits=2)),
    str_cor_plot = paste0("Spearman: R=", round(spearman_rho, 2), 
                          ", p=", format(p_spearman, digits=2))
  )

write_lines(c("", "Part 1: Spearman Correlation with Bootstrapped 95% CI (High vs Low Input)"), outputfile, append=TRUE)
write.table(cor_stats %>% select(-str_cor_plot) %>% mutate(across(starts_with("p_"), format_p_csv)), 
            file=outputfile, append=TRUE, sep=",", row.names=FALSE)


# ==============================================================================
# STATS PART 2: EFFECT OF INPUT AMOUNT (Continuous Linear Model)
# ==============================================================================

anova_stats <- raw_long %>%
  group_by(stat) %>%
  group_modify(~ {
    
    d_mod <- .x %>% 
      select(Plasma_ID, ULP_input_DNA, value) %>% 
      filter(!is.na(value) & !is.na(ULP_input_DNA)) %>%
      mutate(Plasma_ID = factor(Plasma_ID))
    
    # 1. Standard Parametric
    p_anova <- tryCatch({
      m <- lm(value ~ ULP_input_DNA + Plasma_ID, data = d_mod)
      res <- car::Anova(m, type=2) 
      res["ULP_input_DNA", "Pr(>F)"]
    }, error=function(e) {
      m <- lm(value ~ ULP_input_DNA + Plasma_ID, data = d_mod)
      anova(m)["ULP_input_DNA", "Pr(>F)"]
    })
    
    # 2. Permutation
    p_perm <- tryCatch({
      as.numeric(pvalue(independence_test(value ~ ULP_input_DNA | Plasma_ID, data=d_mod, distribution="asymptotic")))
    }, error=function(e) NA)
    
    # 3. Robust (Regression on Differences)
    p_robust <- tryCatch({
      d_diff <- d_mod %>%
        arrange(Plasma_ID, ULP_input_DNA) %>%
        group_by(Plasma_ID) %>%
        summarize(
          d_input = ULP_input_DNA[2] - ULP_input_DNA[1],
          d_val = value[2] - value[1],
          .groups = "drop"
        ) %>%
        filter(!is.na(d_input) & !is.na(d_val))
      
      if (var(d_diff$d_input) < 1e-6) {
        WRS2::yuend(d_diff$d_val, rep(0, nrow(d_diff)), tr=0.2)$p.value
      } else {
        m_rob <- lmrob(d_val ~ d_input, data = d_diff, setting="KS2014")
        summary(m_rob)$coefficients["d_input", "Pr(>|t|)"]
      }
    }, error=function(e) NA)
    
    tibble(
      n_observations = nrow(d_mod),
      p_anova_continuous = p_anova,
      p_perm_continuous = p_perm,
      p_robust_continuous = p_robust
    )
  }) %>%
  mutate(str_eff = paste("Effect of Input (ANOVA): p=", format(p_anova_continuous, digits=2), sep=""))

write_lines(c("", "Part 2: Effect of Continuous Input Amount (Linear Model: y ~ Input + ID)"), outputfile, append=TRUE)
write.table(anova_stats %>% select(-str_eff) %>% mutate(across(starts_with("p_"), format_p_csv)), 
            file=outputfile, append=TRUE, sep=",", row.names=FALSE)


# ==============================================================================
# PLOTTING
# ==============================================================================

theme_replicates <- function(){ 
  theme_cowplot(12) %+replace%
    theme(
      plot.title = element_text(hjust=0,size=8,face="bold",margin=margin(0,0,4,0)),
      plot.subtitle = element_text(hjust=0,size=6,margin=margin(0,0,4,0),lineheight = 1.2), 
      axis.title = element_text(size=7),
      axis.text = element_text(size=6,margin=margin(2,0,0,4)),
      legend.position=c(0.98,0.02),
      legend.justification = c("right", "bottom"),
      legend.margin = margin(3,3,3,3),
      legend.text = element_text(size=4.5,vjust=0.5),legend.key.size = unit(0.20, 'cm'),
      legend.title = element_blank(),
      legend.box.background = element_rect(color="#d9d9d9", size=0.5))
}

# Helper for Labels (Uses str_cor_plot which has NO CI)
get_subtitle <- function(metric) {
  c_str <- cor_stats %>% filter(stat == metric) %>% pull(str_cor_plot)
  e_str <- anova_stats %>% filter(stat == metric) %>% pull(str_eff)
  return(paste(c_str, e_str, sep="\n"))
}

ids <- pairs %>% select(ID) %>% distinct() %>% arrange(ID) %>% pull(ID)

get_shapes <- function(n) {
  base_shapes <- c(0:6, 15:20) 
  if (n <= length(base_shapes)) {
    return(base_shapes[1:n])
  } else {
    return(rep(base_shapes, length.out = n))
  }
}
shape_values <- get_shapes(length(ids))


# --- Plot 1: Tumor Fraction ---
statIn <- "TF"
statNameIn <- "Tumor fraction"
pd <- pairs %>% filter(stat==statIn) 
pd$ID <- factor(pd$ID, levels=ids)

p <- ggplot(pd, aes(x=lower, y=higher))
p <- p + ggtitle(statNameIn, subtitle=get_subtitle(statIn))
p <- p + geom_abline(color="#525252", alpha=0.5)
p <- p + geom_point(aes(shape=ID), color="#cb181d", alpha=0.75, size=1.5)
p <- p + scale_shape_manual(values=shape_values)
p <- p + scale_x_continuous(paste(statNameIn,"(lower amount of DNA)")) 
p <- p + scale_y_continuous(paste(statNameIn,"(higher amount of DNA)")) 
p <- p + theme_replicates() + theme(legend.position="none")
p_TF <- p

# --- Plot 2: Fragment Size ---
statIn <- "Fragment_size_ratio"
statNameIn <- "Fragment size ratio"
pd <- pairs %>% filter(stat==statIn) 
pd$ID <- factor(pd$ID, levels=ids)

max_val <- max(pd$lower, pd$higher, na.rm=T)
min_val <- min(pd$lower, pd$higher, na.rm=T)

p <- ggplot(pd, aes(x=lower, y=higher))
p <- p + ggtitle(statNameIn, subtitle=get_subtitle(statIn))
p <- p + geom_abline(color="#525252", alpha=0.5)
p <- p + geom_point(aes(shape=ID), color="#cb181d", alpha=0.75, size=1.5)
p <- p + scale_shape_manual(values=shape_values)
p <- p + scale_x_continuous(paste(statNameIn,"(lower amount of DNA)"), limits=c(min_val*0.9, max_val*1.35))
p <- p + scale_y_continuous(paste(statNameIn,"(higher amount of DNA)"), limits=c(min_val*0.9, max_val))
p <- p + theme_replicates() 
p_frag <- p

# Assemble
pdfname <- paste0(pdfdir, "Fig_S5.pdf")
grid <- plot_grid(p_TF, p_frag, labels = LETTERS[1:2], nrow=1, label_size = 12, rel_widths = c(1, 1.35))
ggsave(grid, filename=pdfname, width=6.5, height=3)
print(paste("Made", pdfname))