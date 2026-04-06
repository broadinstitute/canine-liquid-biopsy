# ==========================================
# Fragment Length Analysis
# Visualizing cfDNA fragment size distributions by Cancer Type
# ==========================================

library(tidyverse)
library(ggplot2)
library(scales)
# This package automatically finds the root folder of the downloaded project.
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ==============================================================================
# 1. SETUP & FILE PATHS
# ==============================================================================

# Define the input files
infile <- here("data", "all_frag_lengths_OS_LSA_healthy_sample.txt")

# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

print("Loading fragment length data...")
data <- read.delim(infile, stringsAsFactors = FALSE)

# ==========================================
# 2. Preprocessing
# ==========================================

# Use data directly as df_clean (assuming all rows are valid first samples)
df_clean <- data

# Ensure Cancer_type is clean and filter for the 3 groups
target_groups <- c("Osteosarcoma", "Lymphoma", "Healthy")

df_clean <- df_clean %>%
  filter(Cancer_type %in% target_groups) %>%
  mutate(Cancer_type = factor(Cancer_type, levels = target_groups))

print("Samples per group:")
print(table(df_clean %>% select(ID, Cancer_type) %>% distinct() %>% pull(Cancer_type)))

# ==========================================
# 3. Aggregation
# ==========================================
# Calculate the MEAN Proportion for each fragment Size within each Cancer Type.
# This represents the "average profile" for that disease.

df_summary <- df_clean %>%
  group_by(Cancer_type, Size) %>%
  summarize(
    Mean_Prop = mean(Prop, na.rm = TRUE),
    SD_Prop = sd(Prop, na.rm = TRUE),
    SE_Prop = sd(Prop, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Identify the peak (maximum mean proportion) for each Cancer_type
peak_labels <- df_summary %>%
  group_by(Cancer_type) %>%
  filter(Mean_Prop == max(Mean_Prop)) %>%
  ungroup() %>%
  mutate(label_text = paste0(Size, " bp"))

print("Peak fragment sizes:")
print(peak_labels)

# ==========================================
# 4. Plotting
# ==========================================

# Define custom colors
# Osteosarcoma (Orange/Red), Lymphoma (Purple/Blue), Healthy (Green/Grey)
group_colors <- c("Osteosarcoma" = "#d95f02", # Orange
                  "Lymphoma"     = "#7570b3", # Purple
                  "Healthy"      = "#1b9e77") # Green

# --- Common Plot Theme ---
common_theme <- theme_bw() +
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95")
  )

# ---------------------------------------------------------
# Plot 1: Full Range Overlay (Extended to 800bp)
# ---------------------------------------------------------
limit_bp <- 800

p1 <- ggplot(df_summary %>% filter(Size <= limit_bp), 
             aes(x = Size, y = Mean_Prop, fill = Cancer_type, color = Cancer_type)) +
  
  # Use geom_area to create the "filled histogram" look with transparency
  geom_area(position = "identity", alpha = 0.4, size = 0.5) +
  
  # Add labels for the peak size of each group
  geom_text(data = peak_labels, aes(label = label_text), 
            nudge_x = 30, nudge_y = 0.001,
            show.legend = FALSE, size = 3.5, fontface = "bold") +
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(breaks = seq(0, limit_bp, 100), expand = c(0, 0)) +
  scale_y_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0.08))) +
  
  labs(
    title = "cfDNA Fragment Size Distribution (Overlay)",
    subtitle = "Mean proportion of fragments by size across disease groups",
    x = "Fragment Size (bp)",
    y = "Proportion of Fragments",
    fill = "Group",
    color = "Group"
  ) +
  common_theme + 
  theme(legend.position = c(0.85, 0.85))

print(p1)

outfile_1 <- paste0(pdfdir, "Fig_S2_Fragment_Size_Distribution_Overlay.pdf")
ggsave(outfile_1, plot = p1, width = 10, height = 6)
print(paste("Saved plot to:", outfile_1))


# ---------------------------------------------------------
# Plot 2: Zoomed Overlay (0-450bp)
# ---------------------------------------------------------
p2 <- ggplot(df_summary %>% filter(Size >= 0 & Size <= 450), 
             aes(x = Size, y = Mean_Prop, fill = Cancer_type, color = Cancer_type)) +
  
  geom_area(position = "identity", alpha = 0.3) +
  geom_line(size = 1) + 
  
  # Nudged up (y) and right (x) approx 1/3 of grid scales
  geom_text(data = peak_labels, aes(label = label_text), 
            nudge_x = 15, nudge_y = 0.001,
            show.legend = FALSE, size = 4, fontface = "bold") +
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(breaks = seq(0, 450, 50), expand = c(0,0)) +
  scale_y_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0.08))) +
  
  labs(
    title = "cfDNA Fragment Size: Short, Nucleosomal & Dinucleosomal (Zoom)",
    subtitle = "Zoom (0-450bp)",
    x = "Fragment Size (bp)",
    y = "Proportion"
  ) +
  common_theme +
  theme(legend.position = "bottom")

print(p2)

outfile_2 <- paste0(pdfdir, "Fig_S2_Fragment_Size_Distribution_Zoom.pdf")
ggsave(outfile_2, plot = p2, width = 8, height = 6)
print(paste("Saved zoom plot to:", outfile_2))


# ---------------------------------------------------------
# Plot 3: Faceted Full Range
# ---------------------------------------------------------
p3 <- ggplot(df_summary %>% filter(Size <= limit_bp), 
             aes(x = Size, y = Mean_Prop, fill = Cancer_type, color = Cancer_type)) +
  
  geom_area(position = "identity", alpha = 0.6, size = 0.5) +
  
  # Facet Wrap by Cancer Type
  facet_wrap(~Cancer_type, ncol = 1, scales = "free_y") +
  
  # Nudged up (y) and right (x) approx 1/3 of grid scales
  geom_text(data = peak_labels, aes(label = label_text), 
            nudge_x = 30, nudge_y = 0.0015,
            show.legend = FALSE, size = 3.5, fontface = "bold") +
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(breaks = seq(0, limit_bp, 100), expand = c(0, 0)) +
  scale_y_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0.15))) +
  
  labs(
    title = "cfDNA Fragment Size Distribution (Faceted)",
    x = "Fragment Size (bp)",
    y = "Proportion of Fragments"
  ) +
  common_theme + 
  theme(legend.position = "none") # Hide legend since facets explain groups

print(p3)

outfile_3 <- paste0(pdfdir, "Fig_S2_Fragment_Size_Distribution_Faceted.pdf")
ggsave(outfile_3, plot = p3, width = 8, height = 10)
print(paste("Saved faceted plot to:", outfile_3))


# ---------------------------------------------------------
# Plot 4: Faceted Zoom (0-450bp)
# ---------------------------------------------------------
p4 <- ggplot(df_summary %>% filter(Size >= 0 & Size <= 450), 
             aes(x = Size, y = Mean_Prop, fill = Cancer_type, color = Cancer_type)) +
  
  geom_area(position = "identity", alpha = 0.5) +
  geom_line(size = 0.8) + 
  
  facet_wrap(~Cancer_type, ncol = 1, scales = "free_y") +
  
  # Nudged up (y) and right (x) approx 1/3 of grid scales
  geom_text(data = peak_labels, aes(label = label_text), 
            nudge_x = 15, nudge_y = 0.0015,
            show.legend = FALSE, size = 3.5, fontface = "bold") +
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(breaks = seq(0, 450, 50), expand = c(0,0)) +
  scale_y_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0.15))) +
  
  labs(
    title = "cfDNA Fragment Size: Short, Nucleosomal & Dinucleosomal (Faceted)",
    x = "Fragment Size (bp)",
    y = "Proportion"
  ) +
  common_theme +
  theme(legend.position = "none")

print(p4)

outfile_4 <- paste0(pdfdir, "Fig_S2_Fragment_Size_Distribution_Zoom_Faceted.pdf")
ggsave(outfile_4, plot = p4, width = 8, height = 10)
print(paste("Saved faceted zoom plot to:", outfile_4))

# ---------------------------------------------------------
# Plot 5: Sub-nucleosomal High-Res Zoom (0-175bp) with Local Peaks
# ---------------------------------------------------------
# Detect local maxima for the sub-nucleosomal range (e.g., 50bp to 170bp)
# We look for points that are higher than their 2 nearest neighbors on either side
sub_nucleo_limit <- 175

sub_nucleo_peaks <- df_summary %>%
  filter(Size <= sub_nucleo_limit, Size >= 50) %>% 
  group_by(Cancer_type) %>%
  arrange(Size) %>%
  mutate(
    # Check if current point is greater than 2 previous and 2 next points
    is_peak = Mean_Prop > lag(Mean_Prop, 1, default = 0) & 
      Mean_Prop > lead(Mean_Prop, 1, default = 0) &
      Mean_Prop > lag(Mean_Prop, 2, default = 0) & 
      Mean_Prop > lead(Mean_Prop, 2, default = 0)
  ) %>%
  filter(is_peak == TRUE) %>%
  ungroup()

p5 <- ggplot(df_summary %>% filter(Size <= sub_nucleo_limit), 
             aes(x = Size, y = Mean_Prop, fill = Cancer_type, color = Cancer_type)) +
  
  # Minimal area fill
  geom_area(position = "identity", alpha = 0.15) +
  # Strong line to show periodicity clearly
  geom_line(size = 1) + 
  
  # Facet to clearly see the small peaks without overlap
  facet_wrap(~Cancer_type, ncol = 1, scales = "free_y") +
  
  # Add points at detected peaks
  geom_point(data = sub_nucleo_peaks, size = 1.5, alpha = 0.8) +
  
  # Annotate the peaks - using angle to fit more in
  # Nudged up and right approx 1/3 of grid scales
  # X grid is 10bp -> nudge 3bp
  # Y scale max ~0.015 -> nudge ~0.001
  geom_text(data = sub_nucleo_peaks, aes(label = Size), 
            nudge_x = 3, nudge_y = 0.001,
            hjust = 0, # Start text at the nudged point (bottom of text)
            size = 3, angle = 90, show.legend = FALSE, color = "black") +
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  
  # Finer x-axis breaks to see the 10bp periodicity
  scale_x_continuous(breaks = seq(0, sub_nucleo_limit, 10), expand = c(0,0)) +
  scale_y_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0.2))) +
  
  labs(
    title = "Sub-nucleosomal Fragment Size Periodicity (Zoom 0-175bp)",
    subtitle = "Annotated local maxima (peaks) typically observed at ~10bp intervals",
    x = "Fragment Size (bp)",
    y = "Proportion"
  ) +
  common_theme +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = "gray90", linetype = "dashed") # Grid lines to help count periodicity
  )

print(p5)

outfile_5 <- paste0(pdfdir, "Fig_S2_Fragment_Size_Distribution_SubNucleosomal.pdf")
ggsave(outfile_5, plot = p5, width = 8, height = 10)
print(paste("Saved sub-nucleosomal plot to:", outfile_5))