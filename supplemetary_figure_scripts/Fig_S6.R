library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(lubridate)
library(boot) # Required for bootstrapping
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

out_csv  <- paste0(pdfdir, "Table.Fig_S6.stats.csv")
.dir_ensure(pdfdir)

write_lines(
  "Statistical results for analyses presented in Fig. S6 (Spearman + Bootstrap CI)",
  out_csv, append = FALSE
)

# ---------- Helpers ----------
.dir_ensure <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE)

# Scientific Notation Formatter for CSV
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

# Formatter for Plot Labels (simpler)
format_p_plot <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ formatC(p, format = "e", digits = 2),
    TRUE ~ formatC(p, format = "f", digits = 3)
  )
}

# Bootstrap Function
boot_spearman <- function(data, indices) {
  d <- data[indices,]
  return(cor(d$x, d$y, method = "spearman", use = "complete.obs"))
}

# Metric labels / ordering
rename_stats <- tibble(
  order       = c(2, 1, 3),
  statName    = c("Tumor fraction", "cfDNA concentration", "Fragment size ratio"),
  axisLabel = c("Tumor fraction", "cfDNA concentration (ng/mL)", "Fragment size ratio"),
  stat      = c("TF", "log_DNA_conc", "Fragment_size_ratio")
) %>% arrange(order)

# ---------- Data Loading & Cleaning ----------
raw <- as_tibble(read.delim(
  infile,
  header = TRUE,
  sep = "\t",
  na.strings = c("NA","N/A","","UNK","?")
))
raw <- raw %>% filter(!is.na(Sample_ID))
if ("X" %in% names(raw)) raw$X <- NULL

# Dates
raw <- raw %>%
  mutate(Date_of_sample = str_replace(Date_of_sample, " \\(approx\\)", "")) %>%
  mutate(Date_of_sample = parse_date_time(Date_of_sample, "m/d/y"))

# Status
raw <- raw %>%
  mutate(Cancer_status = if_else(Cancer_type != "Healthy", "cancer", "healthy"))

# Log DNA concentration (non-positive → NA)
raw <- raw %>%
  mutate(log_DNA_conc = if_else(
    is.finite(DNA_conc) & DNA_conc > 0,
    log10(DNA_conc),
    NA_real_
  ))

# Keep first samples, one rep
long <- raw %>%
  filter(First_Sample, is.na(replicate) | replicate == "Rep1") %>%
  select(ID, Sample_ID, TF, Fragment_size_ratio, log_DNA_conc, Cancer_status, Cancer_type, Plasma_volume_estimated) %>%
  pivot_longer(
    c(log_DNA_conc, TF, Fragment_size_ratio),
    names_to = "stat",
    values_to = "value"
  ) %>%
  filter(is.finite(value)) %>%
  left_join(rename_stats %>% select(stat, statName, order), by = "stat")

# Build metric pairs per sample (statX, statY) with orderY > orderX
pairs <- long %>%
  select(ID, Sample_ID, Cancer_status, Cancer_type, Plasma_volume_estimated,
         statX = stat, valueX = value, statNameX = statName, orderX = order) %>%
  inner_join(
    long %>%
      select(ID, Sample_ID, Cancer_status, Cancer_type, Plasma_volume_estimated,
             statY = stat, valueY = value, statNameY = statName, orderY = order),
    by = c("ID","Sample_ID","Cancer_status","Cancer_type", "Plasma_volume_estimated"),
    relationship = "many-to-many"
  ) %>%
  filter(orderY > orderX) %>%
  drop_na(valueX, valueY)

# Grouping: (a) status (healthy vs cancer), (b) cancer types (lymphoma vs other cancers)
groups <- raw %>%
  transmute(ID, group = "status", status = Cancer_status) %>%
  distinct() %>%
  bind_rows(
    raw %>%
      filter(Cancer_status == "cancer") %>%
      transmute(
        ID,
        group  = "types",
        status = if_else(Cancer_type == "Lymphoma", "lymphoma", "other cancers")
      ) %>%
      distinct()
  )

pairs <- pairs %>%
  left_join(groups, by = "ID", relationship = "many-to-many")

# Filter out the generic "cancer" group for plotting (we split by type instead)
pdAll <- pairs %>%
  mutate(CaseGroup = status) %>%
  filter(CaseGroup != "cancer")

# Using fillable shapes: 21 (Circle), 25 (Inv Triangle)
col_map   <- c("healthy" = "#525252", "lymphoma" = "#ef3b2c", "other cancers" = "#8073ac")
shape_map <- c("healthy" = 21, "lymphoma" = 21, "other cancers" = 25)

# Special map for Panel A (includes "cancer" aggregate)
col_map_A <- c("healthy" = "#525252", "cancer" = "#a50f15", 
               "lymphoma" = "#ef3b2c", "other cancers" = "#8073ac")

# ==============================================================================
# STATISTICS: SPEARMAN + BOOTSTRAP CI
# ==============================================================================

cor_diag <- pairs %>%
  # Group by comparison type and specific pair
  group_by(group, status, statNameX, statNameY, statX, statY, orderX, orderY) %>%
  group_modify(~{
    # Prepare vector for correlation
    df_cor <- tibble(x = .x$valueX, y = .x$valueY) %>% drop_na()
    n <- nrow(df_cor)
    
    if (n < 3) {
      return(tibble(N=n, rho=NA, p=NA, ci_low=NA, ci_high=NA))
    }
    
    # 1. Standard Spearman
    ct <- cor.test(df_cor$x, df_cor$y, method = "spearman", exact = FALSE)
    
    # 2. Bootstrap CI
    set.seed(42)
    boot_res <- boot(data = df_cor, statistic = boot_spearman, R = 2000)
    # Handle errors in boot.ci (e.g. if all values are identical)
    ci <- tryCatch({
      boot.ci(boot_res, type = "perc")$percent[4:5]
    }, error = function(e) c(NA, NA))
    
    tibble(
      N = n,
      rho = ct$estimate,
      p = ct$p.value,
      ci_low = ci[1],
      ci_high = ci[2]
    )
  }) %>%
  ungroup()

# Write to CSV
csv_out <- cor_diag %>%
  mutate(across(c(p, rho, ci_low, ci_high), ~format_p_csv(.))) %>%
  select(group, status, Pair1=statNameX, Pair2=statNameY, N, rho, ci_low, ci_high, p)

# Suppress warnings for column names append
suppressWarnings(
  write.table(csv_out, file = out_csv, append = TRUE, sep = ",", row.names = FALSE)
)


# ==============================================================================
# PLOTTING: SCATTER PANELS (Row 1)
# ==============================================================================

make_scatter <- function(df_all, xvar, yvar, stats_df) {
  pd <- df_all %>% filter(statX == xvar, statY == yvar)
  if (nrow(pd) == 0) return(ggplot() + theme_void())
  
  # Prepare flag for Estimated vs Not Estimated
  # Logic: Apply to ALL plots as requested (removed xvar condition)
  pd <- pd %>% 
    mutate(
      is_est = (isTRUE(Plasma_volume_estimated) | Plasma_volume_estimated == "TRUE"),
      VolStatus = if_else(is_est, "Plasma estimated", "Not estimated")
    )
  
  # Ensure factors for consistent legend ordering
  pd$VolStatus <- factor(pd$VolStatus, levels = c("Plasma estimated", "Not estimated"))
  
  # Axis labels
  xlabel <- rename_stats %>% filter(stat == xvar) %>% pull(axisLabel) %>% unique()
  ylabel <- rename_stats %>% filter(stat == yvar) %>% pull(axisLabel) %>% unique()
  title  <- paste0(
    (rename_stats %>% filter(stat == xvar) %>% pull(statName) %>% unique()),
    " vs.\n",
    tolower(rename_stats %>% filter(stat == yvar) %>% pull(statName) %>% unique())
  )
  
  # y-positioning for labels
  y_min <- min(pd$valueY, na.rm = TRUE)
  y_max <- max(pd$valueY, na.rm = TRUE)
  y_inc <- (y_max - y_min) / 15
  
  # Order labels consistently
  levs <- intersect(c("healthy","lymphoma","other cancers"), unique(pd$CaseGroup))
  
  # Create label dataframe from the pre-calculated stats
  lab_df <- tibble(CaseGroup = levs) %>%
    left_join(stats_df %>% filter(statX == xvar, statY == yvar), by = c("CaseGroup" = "status")) %>%
    mutate(
      lab = paste0("rho = ", round(rho, 3), "; p = ", format_p_plot(p)),
      x_pos = min(pd$valueX, na.rm = TRUE),
      y_pos = y_max + (length(levs) - row_number() + 1) * y_inc
    )
  
  # Determine Y-Scale ONCE
  if (yvar == "TF") {
    y_scale <- scale_y_continuous(
      ylabel,
      expand = expansion(mult = c(0.05, 0.15)),
      breaks = c(0, 0.5, 1)
    )
  } else {
    y_scale <- scale_y_continuous(
      ylabel,
      expand = expansion(mult = c(0.05, 0.15))
    )
  }
  
  # Base plot
  p <- ggplot(pd, aes(x = valueX, y = valueY)) +
    
    # Layer 1: Estimated Plasma (Filled with Group Color)
    # Mapping alpha to VolStatus creates the legend key for "Plasma estimated"
    geom_point(data = pd %>% filter(is_est),
               aes(color = CaseGroup, shape = CaseGroup, fill = CaseGroup, alpha = VolStatus),
               size = 1.2) +
    
    # Layer 2: Not Estimated Plasma (Filled with NA = Hollow)
    # Mapping alpha to VolStatus creates the legend key for "Not estimated"
    geom_point(data = pd %>% filter(!is_est),
               aes(color = CaseGroup, shape = CaseGroup, alpha = VolStatus),
               fill = NA, 
               size = 1.2) +
    
    # Regression lines (Hidden from legend)
    geom_smooth(
      aes(group = CaseGroup, color = CaseGroup, fill = CaseGroup),
      method = "lm", se = TRUE,
      alpha = 0.15, linewidth = 0.5,
      show.legend = FALSE
    ) +
    
    # Text (Hidden from legend)
    geom_text(
      data = lab_df,
      aes(x = x_pos, y = y_pos, label = lab, color = CaseGroup),
      hjust = 0, vjust = 0, size = 2.2,
      show.legend = FALSE
    ) +
    
    # Scales
    scale_color_manual(values = col_map, breaks = names(col_map)) +
    scale_fill_manual(values  = col_map, breaks = names(col_map)) +
    scale_shape_manual(values = shape_map, breaks = names(shape_map)) +
    
    # Alpha scale used for Legend Annotation of Volume Status
    # drop=FALSE to prevent error when data doesn't have both levels
    scale_alpha_manual(values = c("Plasma estimated" = 0.6, "Not estimated" = 0.6), 
                       name = NULL,
                       drop = FALSE) +
    
    y_scale +
    
    # Unified Legend Guides
    guides(
      # Consolidated main overrides into 'color' to fix warnings
      color = guide_legend(override.aes = list(alpha = 1, size = 2.5), order = 1),
      fill  = guide_legend(order = 1),
      shape = guide_legend(order = 1),
      
      # Create visual distinction in legend for Filled vs Hollow
      # Since drop=FALSE, we are guaranteed 2 keys, so this 2-value vector won't error
      alpha = guide_legend(override.aes = list(shape = 21, fill = c("black", NA), color = "black"), order = 2)
    ) +
    
    ggtitle(title) +
    theme_cowplot(12) +
    theme(
      plot.title    = element_text(hjust = 0, size = 8, face = "bold",
                                   margin = margin(0,0,4,0)),
      plot.subtitle = element_blank(),
      legend.position = "bottom",
      legend.text   = element_text(size = 6),
      legend.title  = element_blank(),
      axis.title    = element_text(size = 7),
      axis.text     = element_text(size = 6),
      legend.key.size = unit(0.22, "cm")
    )
  
  # Special x scales by metric
  if (xvar == "log_DNA_conc") {
    p <- p + scale_x_continuous(
      xlabel,
      breaks = -1:3,
      labels = 10^(-1:3)
    )
  } else {
    p <- p + scale_x_continuous(xlabel)
  }
  
  p
}

# Build three scatter panels
p_cor12 <- make_scatter(pdAll, xvar = "log_DNA_conc",   yvar = "TF", cor_diag)
p_cor13 <- make_scatter(pdAll, xvar = "log_DNA_conc",   yvar = "Fragment_size_ratio", cor_diag)
p_cor23 <- make_scatter(pdAll, xvar = "TF",             yvar = "Fragment_size_ratio", cor_diag)

# EXTRACT SHARED LEGEND from the first plot (Panel B)
# This legend contains both Cancer Status and Plasma Volume keys
shared_legend <- get_legend(p_cor12 + theme(legend.box.margin = margin(0, 0, 0, 12)))

# Remove legends from individual plots to avoid duplication
p_cor12 <- p_cor12 + theme(legend.position = "none")
p_cor13 <- p_cor13 + theme(legend.position = "none")
p_cor23 <- p_cor23 + theme(legend.position = "none")

# Combine plots into row
row_cors <- plot_grid(
  p_cor12, p_cor13, p_cor23,
  labels = LETTERS[2:4], ncol = 3, label_size = 12
)

# Attach shared legend to the bottom of the plots
row_cors_w_legend <- plot_grid(row_cors, shared_legend, ncol = 1, rel_heights = c(1, 0.1))

# ==============================================================================
# PLOTTING: FOREST PLOT SUMMARY (Row 2/Panel A)
# ==============================================================================

# Build facets for Forest Plot
facets <- cor_diag %>%
  distinct(statNameX, statNameY, orderX, orderY) %>%
  mutate(
    facet       = paste0(statNameX, " vs.\n", tolower(statNameY)),
    facet_order = orderX + orderY/10
  )

# Prepare summary data
cor_summary <- cor_diag %>%
  left_join(facets, by = c("statNameX","statNameY","orderX","orderY")) %>%
  mutate(
    pstr  = paste0("p = ", format_p_plot(p)),
    facet = factor(
      facet,
      levels = facets %>% arrange(facet_order) %>% pull(facet)
    )
  )

# y placement helper
ylabs <- cor_summary %>%
  distinct(group, status) %>%
  arrange(group, status) %>%
  mutate(ypos = ave(as.integer(factor(status)), group, FUN = seq_along))

# Plot B (Forest Plot style)
pB <- cor_summary %>%
  left_join(ylabs, by = c("group","status")) %>%
  ggplot(aes(y = status, x = rho)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.25, color = "#525252") +
  # Error Bars (Bootstrap CI)
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = status, color = status)
  ) +
  # Point Estimate (Rho)
  geom_point(
    aes(color = status, fill = status),
    shape = 16,
    data = ~ dplyr::filter(.x, !is.na(p) & is.finite(rho)),
    size = 1
  ) +
  # Empty circle for missing data
  geom_point(
    color = "#FFFFFF", shape = 16,
    data = ~ dplyr::filter(.x, is.na(p) | !is.finite(rho)),
    size = 2.5
  ) +
  geom_point(
    aes(color = status),
    fill = NA, shape = 21, size = 1
  ) +
  # P-value Text
  geom_text(
    aes(x = rho, y = ypos + 0.5, label = pstr, color = status),
    hjust = 0.5, vjust = 1, size = 1.85
  ) +
  scale_x_continuous("correlation (Spearman rho)", expand = c(0, 0.5)) +
  scale_y_discrete("", expand = expansion(add = 0.3)) +
  facet_grid(
    rows = vars(group),
    cols = vars(facet),
    scales = "free",
    space  = "free_y"
  ) +
  # Use the explicit named color map (col_map_A)
  scale_color_manual(values = col_map_A) +
  scale_fill_manual(values  = col_map_A) +
  ggtitle("Correlation (Spearman) between liquid biopsy metrics in healthy and cancer dogs") +
  theme_cowplot(12) +
  theme(
    plot.title        = element_text(hjust = 0, size = 8, face = "bold",
                                     margin = margin(0,0,4,0)),
    legend.position   = "none",
    axis.title.y      = element_blank(),
    axis.title.x      = element_text(size = 7),
    axis.text.y       = element_text(size = 6, margin = margin(0,0,0,4)),
    axis.text.x       = element_text(size = 7),
    strip.text.x      = element_text(size = 6, hjust = 0.5,
                                     face = "bold", color = "#f7f7f7"),
    strip.background.x = element_rect(fill = "#525252"),
    strip.text.y      = element_blank(),
    strip.background.y = element_blank()
  )

row_cor_range <- plot_grid(pB, labels = LETTERS[1], ncol = 1, label_size = 12)

# ---------- Export figure ----------
pdfname <- paste0(pdfdir, "Fig_S6.pdf")
grid <- plot_grid(
  row_cor_range,
  row_cors_w_legend,
  labels = NULL,
  ncol = 1,
  rel_heights = c(0.65, 1)
)
ggsave(grid, filename = pdfname, width = 6.5, height = 4.5)
message(paste("Made", pdfname))