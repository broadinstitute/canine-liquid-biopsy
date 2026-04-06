# =========================
# STUBES checks (integrated with Robust Permutation Tests)
# =========================
library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(corrr)
library(lubridate)
library(coin) 
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

outputfile <- paste0(pdfdir,"Table.Fig_S14.stats.csv")
write_lines("Statistical results for analyses presented in Fig. S14", outputfile, append = FALSE)

# =========================
# Helpers & Formatters
# =========================

# CHANGED: Exact P-value formatting for Plot
format_p_plot <- function(p) {
  sapply(p, function(val) {
    if (is.na(val)) return("NA")
    
    # Handle absolute machine zero (rare, but possible with asymptotic tests)
    if (val < 2.2e-16) return("< 2.2e-16")
    
    # If small, use scientific notation (e.g. 1.2e-05)
    if (val < 0.001) {
      formatC(val, format = "e", digits = 2)
    } else {
      # Otherwise standard decimal (e.g. 0.042)
      formatC(val, format = "f", digits = 3)
    }
  })
}

# Helper to extract p-value from coin object safely
# Explicit coin::pvalue call to avoid namespace conflicts
get_coin_p <- function(x) {
  tryCatch({
    as.numeric(coin::pvalue(x))
  }, error = function(e) NA)
}

# 1. Independent-groups Permutation (Wilcoxon rank-sum)
.perm_p_wilcox_indep <- function(x, g, B = 10000, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = x, g = as.factor(g))
  
  # Wrapped in tryCatch and explicit coin:: calls
  tryCatch({
    test <- coin::wilcox_test(x ~ g, data = d, distribution = coin::approximate(nresample = B))
    get_coin_p(test)
  }, error = function(e) {
    NA
  })
}

# 2. Paired Permutation (Signed-Rank)
.perm_p_wilcox_paired <- function(x1, x2, B = 10000, seed = 1) {
  set.seed(seed)
  ok <- is.finite(x1) & is.finite(x2)
  x1 <- x1[ok]; x2 <- x2[ok]
  
  if(length(x1) < 2) return(NA)
  
  d <- data.frame(
    val = c(x1, x2),
    group = factor(rep(c("A", "B"), each = length(x1))),
    id = factor(rep(1:length(x1), 2))
  )
  
  # Wrapped in tryCatch and explicit coin:: calls
  tryCatch({
    test <- coin::wilcoxsign_test(val ~ group | id, data = d, distribution = coin::approximate(nresample = B))
    get_coin_p(test)
  }, error = function(e) {
    NA
  })
}

# Brown–Forsythe (median-centered Levene)
.brown_forsythe <- function(x, g) {
  g <- as.factor(g)
  z <- abs(x - ave(x, g, FUN = function(v) median(v, na.rm = TRUE)))
  fit <- lm(z ~ g); a <- anova(fit)
  list(p = a$`Pr(>F)`[1])
}

# Diagnostics for unpaired comparisons
check_t_independent <- function(x, g, alpha = 0.05, B_perm = 10000, seed = 1) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(as.factor(g[ok]))
  if (nlevels(g) != 2) stop("Group factor must have exactly 2 levels.")
  
  tt_welch <- t.test(x ~ g, var.equal = FALSE)
  wt_std <- wilcox.test(x ~ g, exact = FALSE)
  p_perm <- .perm_p_wilcox_indep(x, g, B = B_perm, seed = seed)
  bf <- .brown_forsythe(x, g)
  
  list(
    welch_p = tt_welch$p.value,
    wilcox_p = wt_std$p.value,
    perm_p = p_perm, 
    levene_bf_p = bf$p
  )
}

# Build manual p-value labels for stat_pvalue_manual()
make_wilcox_labels <- function(df, x = "xset", y = "value", pairs,
                               step_mult = 0.06) {
  if (!all(c(x, y) %in% names(df))) return(tibble())
  pairs <- purrr::keep(pairs, ~ all(.x %in% unique(df[[x]])))
  if (length(pairs) == 0) return(tibble())
  
  y_max <- max(df[[y]], na.rm = TRUE)
  y_min <- min(df[[y]], na.rm = TRUE)
  y_span <- y_max - y_min
  ymax_by_group <- df %>% group_by(.data[[x]]) %>% summarise(ymax = max(.data[[y]], na.rm = TRUE), .groups = "drop")
  
  purrr::map2_dfr(pairs, seq_along(pairs), function(pr, i) {
    dsub <- df %>% dplyr::filter(.data[[x]] %in% pr) %>% droplevels()
    
    # Calculate Standard Wilcoxon P-value
    f <- as.formula(paste(y, "~", x))
    p <- wilcox.test(f, data = dsub, exact = FALSE)$p.value
    
    y_base <- max(ymax_by_group$ymax[match(pr, ymax_by_group[[x]])], na.rm = TRUE)
    tibble(
      group1 = pr[1],
      group2 = pr[2],
      # Apply formatter (Exact P)
      p.label = paste0("p=", format_p_plot(p)), 
      y.position = y_base + i * step_mult * y_span
    )
  })
}

# =========================
# Analysis
# =========================

rename_stats <- tibble(
  order     = c(2,1,3),
  title     = c("Tumor fraction","cfDNA concentration","Fragment size ratio"),
  statistic = c("Tumor fraction","cfDNA concentration (ng/mL)","Fragment size ratio"),
  statName  = c("TF","log_DNA_conc","Fragment_size_ratio")
) %>% arrange(order)

# Tab-delimited input
raw <- as_tibble(read.delim(infile, header = TRUE, sep = "\t", na.strings = c("NA","N/A","","UNK")))
raw <- raw %>% mutate(Sample_type = str_remove(Sample_type, "_plasma")) %>% filter(!is.na(Sample_type))

# add lymphoma flag
raw <- raw %>% mutate(Lymphoma = Cancer_type == "Lymphoma")

# date parse
raw <- raw %>%
  mutate(Date_of_sample = str_replace(Date_of_sample, " \\(approx\\)", "")) %>%
  mutate(Date_of_sample = parse_date_time(Date_of_sample, "m/d/y"))

# recode healthy/cancer
raw <- raw %>% mutate(Cancer_status = if_else(Cancer_type != "Healthy", "cancer", "healthy"))

# Keep cancer samples only
raw <- raw %>% filter(Cancer_status == "cancer")

# keep + long pivot
raw <- raw %>%
  select(ID, First_Sample, Sample_type, Cancer_status, Lymphoma, Date_of_sample, DNA_conc, TF, Fragment_size_ratio) %>%
  mutate(log_DNA_conc = log10(DNA_conc)) %>%
  select(-DNA_conc) %>%
  pivot_longer(c(log_DNA_conc, TF, Fragment_size_ratio), names_to = "statName", values_to = "value") %>%
  filter(!is.na(value))

# ------------------------------
# PAIRED: robust detection
# ------------------------------
paired_keys <- raw %>%
  distinct(ID, Date_of_sample, statName, Sample_type) %>%
  count(ID, Date_of_sample, statName) %>%
  filter(n >= 2) %>%
  select(ID, Date_of_sample, statName)

reps <- raw %>%
  semi_join(paired_keys, by = c("ID","Date_of_sample","statName")) %>%
  select(ID, Date_of_sample, statName, Sample_type, value)

reps_wide <- reps %>%
  pivot_wider(names_from = Sample_type, values_from = value) %>%
  filter(is.finite(EDTA), is.finite(Streck))

# Paired Stats
paired <- reps_wide %>%
  group_by(statName) %>%
  summarise(
    n_pairs     = n(),
    mean_EDTA   = mean(EDTA, na.rm = TRUE),
    mean_Streck = mean(Streck, na.rm = TRUE),
    p_ttest     = t.test(EDTA, Streck, paired = TRUE)$p.value,
    p_wilcox    = wilcox.test(EDTA, Streck, paired = TRUE)$p.value,
    p_perm      = .perm_p_wilcox_paired(EDTA, Streck),
    .groups     = "drop"
  ) %>%
  # UPDATED: Add FDR Adjusted P-value
  mutate(
    p_wilcox_fdr = p.adjust(p_wilcox, method = "fdr"),
    # Label for plot using Standard Paired Wilcoxon (Exact P)
    subtitle = paste0(
      "Paired Wilcoxon p=", format_p_plot(p_wilcox),      
      "\nPaired means: EDTA = ", round(mean_EDTA, 2),
      " & Streck = ", round(mean_Streck, 2)
    )
  )

# ------------------------------
# UNPAIRED pool
# ------------------------------
paired_ids <- unique(paired_keys$ID)

unpaired <- raw %>%
  filter(First_Sample, !(ID %in% paired_ids)) %>%
  select(-any_of(c("First_Sample","Cancer_status","Date_of_sample"))) %>%
  ungroup()

# sets
sets  <- unpaired %>% select(ID) %>% distinct() %>% mutate(set = "all")
sets <- unpaired %>% filter(!Lymphoma) %>% select(ID) %>% distinct() %>% mutate(set = "no lymphomas") %>% bind_rows(sets)
sets <- unpaired %>% filter(Lymphoma) %>% select(ID) %>% distinct() %>% mutate(set = "lymphomas") %>% bind_rows(sets)
sets <- reps %>% ungroup() %>% select(ID) %>% distinct() %>% mutate(set = "reps") %>% bind_rows(sets)

# ===========================
# Summary tables to CSV
# ===========================

# B. Basic summaries
table <- unpaired %>%
  left_join(sets, by = "ID") %>%
  group_by(Sample_type, set, statName) %>%
  summarise(
    N = n(),
    min = min(value), max = max(value),
    median = median(value), mean = mean(value), sd = sd(value),
    .groups = "drop"
  )
write_lines(c("", "B. Comparison of liquid biopsy values by tube type"), outputfile, append = TRUE)
write.table(table, file = outputfile, append = TRUE, sep = ",", row.names = FALSE)

# C. Unpaired diagnostics
diag_unpaired <- unpaired %>%
  left_join(sets, by = "ID") %>%
  group_by(statName, set) %>%
  group_modify(~{
    d <- .x %>% filter(Sample_type %in% c("EDTA","Streck")) %>% droplevels()
    if (n_distinct(d$Sample_type) < 2) return(tibble())
    
    ck <- check_t_independent(d$value, d$Sample_type, B_perm = 10000, seed = 1)
    
    tibble(
      n_EDTA        = sum(d$Sample_type == "EDTA"),
      n_Streck      = sum(d$Sample_type == "Streck"),
      p_welch       = ck$welch_p,
      p_wilcox      = ck$wilcox_p,
      p_perm_wilcox = ck$perm_p,
      p_levene_bf   = ck$levene_bf_p
    )
  }) %>% 
  ungroup() %>%
  # Add FDR Adjusted P-value
  mutate(p_wilcox_fdr = p.adjust(p_wilcox, method = "fdr"))

write_lines(c("", "C. Unpaired Diagnostics: Welch vs Wilcoxon vs Permutation"), outputfile, append = TRUE)
write.table(diag_unpaired, file = outputfile, append = TRUE, sep = ",", row.names = FALSE)

# D. Lymphoma vs Non-Lymphoma within Streck
diag_streck_lymph <- unpaired %>%
  dplyr::filter(Sample_type == "Streck") %>%
  dplyr::group_by(statName) %>%
  dplyr::group_modify(~{
    d <- .x %>% dplyr::filter(!is.na(Lymphoma))
    if (dplyr::n_distinct(d$Lymphoma) < 2) return(tibble())
    g <- factor(ifelse(d$Lymphoma, "Lymphoma", "Non-Lymphoma"), levels = c("Non-Lymphoma","Lymphoma"))
    
    ck <- check_t_independent(d$value, g, B_perm = 10000, seed = 1)
    
    tibble(
      n_non_lymph   = sum(!d$Lymphoma),
      n_lymph       = sum( d$Lymphoma),
      p_welch       = ck$welch_p,
      p_wilcox      = ck$wilcox_p,
      p_perm_wilcox = ck$perm_p,
      p_levene_bf   = ck$levene_bf_p
    )
  }) %>% 
  dplyr::ungroup() %>%
  # Add FDR Adjusted P-value
  mutate(p_wilcox_fdr = p.adjust(p_wilcox, method = "fdr"))

write_lines(c("", "D. Lymphoma vs Non-Lymphoma within Streck (Unpaired Diagnostics)"), outputfile, append = TRUE)
write.table(diag_streck_lymph, file = outputfile, append = TRUE, sep = ",", row.names = FALSE)

# E. Paired Summary
write_lines(c("", "E. Paired comparisons (T-test vs Wilcoxon vs Permutation)"), outputfile, append = TRUE)
write.table(paired, file = outputfile, append = TRUE, sep = ",", row.names = FALSE)

# ===========================
# Plotting
# ===========================

theme_compare <- function(){ 
  theme_cowplot(12) %+replace%
    theme(
      plot.title = element_text(hjust = 0, size = 8, face = "bold", margin = margin(0,0,4,0)),
      plot.subtitle = element_text(hjust = 0, size = 6, margin = margin(0,0,4,0), lineheight = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8, angle = 90),
      axis.text.x  = element_text(size = 5),
      axis.text.y  = element_text(size = 6),
      legend.position = c(0.98,0.02),
      legend.justification = c("right","bottom"),
      legend.text = element_text(size = 4.5, vjust = 0.5),
      legend.key.size = unit(0.18, 'cm'),
      legend.title = element_blank()
    )
}

pdAll <- unpaired %>%
  select(-any_of("Lymphoma")) %>% mutate(subset = "all") %>%
  left_join(sets, by = "ID", relationship = "many-to-many") %>%
  bind_rows(reps %>% mutate(set = "all", subset = "paired"))
pdAll <- pdAll %>% mutate(xset = paste(set, Sample_type))
pdAll <- pdAll %>%
  select(xset, ID) %>% distinct() %>% group_by(xset) %>% count() %>%
  filter(n >= 10) %>% select(-n) %>% inner_join(pdAll, by = "xset")
pdAll <- pdAll %>% left_join(rename_stats, by = "statName")
pdAll$set <- factor(pdAll$set, levels = c("all","no lymphomas","lymphomas"))

pairs_for_plots <- list(
  c("all EDTA",             "all Streck"),
  c("no lymphomas Streck", "no lymphomas EDTA"),
  c("no lymphomas Streck", "lymphomas Streck")
)

# ---- Plot Generator ----
make_stub_plot <- function(metric) {
  pd <- pdAll %>% filter(statName == metric & subset == "all")
  
  # Labels
  xlabels <- pd %>% group_by(Sample_type, set, xset) %>% count() %>%
    mutate(xlabel = paste0(Sample_type, "\n(N=", n, ")") %>% str_replace("\\_", "\n"))
  xlabels <- xlabels %>% left_join(tibble(set = c("all","no lymphomas","lymphomas"), order = c(1,2,3)), by = "set") %>%
    mutate(order = if_else(Sample_type == "Streck", order + 0.5, order))
  pd$xset <- factor(pd$xset, levels = xlabels %>% arrange(order) %>% pull(xset))
  
  title <- paste0("Streck vs EDTA\n(", unique(pd$title), ")")
  
  # Subtitle from Paired Stats
  sub_text <- paired %>% filter(statName == metric) %>% pull(subtitle)
  subtitle <- paste0("Dogs with cancer\n", sub_text)
  
  # Use Standard Wilcoxon for labels with Exact Formatting
  perm_lbls <- make_wilcox_labels(pd, x = "xset", y = "value", pairs = pairs_for_plots, step_mult = 0.08)
  
  # Paired points background
  pdReps <- pdAll %>% filter(statName == metric & subset == "paired")
  lines <- pdReps %>% filter(str_detect(Sample_type, "EDTA")) %>% select(ID, xset, value) %>%
    full_join(pdReps %>% filter(str_detect(Sample_type, "Streck")) %>% select(ID, xset, value) %>% rename(xend = xset, yend = value), by = "ID")
  
  p <- ggplot(pd, aes(x = xset, y = value)) +
    geom_boxplot(aes(color = set), fill = "#FFFFFF", alpha = 0.75, outlier.size = 0.75, outlier.shape = NA) +
    ggpubr::stat_pvalue_manual(perm_lbls, label = "p.label", xmin = "group1", xmax = "group2",
                               y.position = "y.position", tip.length = 0, size = 2.3, bracket.size = 0.3) +
    geom_jitter(aes(color = set), size = 0.75, alpha = 0.5, height = 0, width = 0.2, shape = 21) +
    geom_point(data = pdReps, color = "#252525", size = 1.5, alpha = 0.75, shape = 16) +
    scale_color_manual(values = c("#a50f15","#8073ac","#ef3b2c")) +
    geom_segment(aes(xend = xend, yend = yend), data = lines, color = "#252525", alpha = 1, linewidth = 0.5) +
    scale_x_discrete("Tube type", breaks = xlabels$xset, labels = xlabels$xlabel) +
    ggtitle(title, subtitle = subtitle) +
    theme_compare()
  
  # Scales
  if(metric == "log_DNA_conc") {
    p <- p + scale_y_continuous(unique(pd$statistic), expand = expansion(mult = c(0.2,0.1)), breaks = c(-1:3), labels = 10**c(-1:3))
  } else if (metric == "TF") {
    p <- p + scale_y_continuous("tumor fraction", expand = expansion(mult = c(0.2,0.1)), breaks = c(0,0.5,1))
  } else {
    p <- p + scale_y_continuous("fragment size ratio", expand = expansion(mult = c(0.2,0.1)))
  }
  
  return(p)
}

p_tube_conc <- make_stub_plot("log_DNA_conc")
p_tube_TF <- make_stub_plot("TF")
p_tube_frag <- make_stub_plot("Fragment_size_ratio")

# ---- Export figure ----
pdfname <- paste0(pdfdir,"Fig_S14.pdf")
grid <- plot_grid(p_tube_conc, p_tube_TF, p_tube_frag, labels = LETTERS[1:3], ncol = 3, label_size = 12)
ggsave(grid, filename = pdfname, width = 6.5, height = 3)
message(paste("Made", pdfname))