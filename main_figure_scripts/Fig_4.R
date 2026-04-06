library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(corrr)
library(car)
library(boot) # Added for CIs
library(lubridate)
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

outputfile <- paste0(pdfdir, "Table.Fig_4.stats.csv")

write_lines(c("Statistical results for analyses presented in Fig. 4"), outputfile, append=FALSE)

rename_stats <- tibble(
  order=c(2,1,3),
  statistic=c("Tumor fraction","cfDNA concentration","Fragment size ratio"),
  stat=c("TF","log_DNA_conc","Fragment_size_ratio") 
) %>% arrange(order)

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
    p < 0.0001 ~ formatC(p, format = "e", digits = 2),
    TRUE ~ formatC(p, format = "f", digits = 4)
  )
}

# --- Data Loading ---
raw <- as_tibble(read.csv(infile,header=T,sep="\t",na.strings = c("NA","N/A","","UNK")))
cbc <- as_tibble(read.csv(infile_cbc,header=T,sep="\t",na.strings = c("NA","N/A","","UNK")))

raw <- raw %>% mutate(Sample_type=str_remove(Sample_type,"_plasma"))

# format Date_of_sample column as a date
raw <- raw %>% mutate(Date_of_sample=str_replace(Date_of_sample," \\(approx\\)","")) %>% mutate(Date_of_sample=parse_date_time(Date_of_sample,"m/d/y"))
raw <- raw %>% mutate(Cancer_status=if_else(Cancer_type!="Healthy","cancer","healthy"))

dogs <-  raw %>% select(ID,Cancer_status,Cancer_type,Sex,Breed) %>% distinct()
raw <- raw %>% select(-Cancer_status,-Cancer_type,-Sex,-Breed) %>% distinct()
# add column for mixed or single breed
dogs <- dogs %>% filter(!is.na(Breed)) %>% mutate(Ancestry=if_else(Breed=="mixed breed","mixed","single")) %>% select(ID,Ancestry) %>% distinct() %>% right_join(dogs)

breeds <- dogs %>% select(ID,Breed) %>% distinct() %>% group_by(Breed) %>% count() %>% mutate(Breed_v2=if_else(n>=2,Breed,"other breeds")) %>% rename(ndogs_in_breed=n) %>% select(Breed,ndogs_in_breed,Breed_v2) %>% distinct()
dogs <- dogs %>% left_join(breeds)

# recode Cancer_status as healthy/cancer
dogs <- dogs %>% mutate(Cancer_status=if_else(Cancer_status=="yes","cancer",if_else(Cancer_status=="no","healthy","other")))

# add column for whether cancer type is X or something else 
dogs <- dogs %>% mutate(Lymphoma=if_else(Cancer_type=="Lymphoma",TRUE,FALSE))

okbreeds_gt1 <- dogs %>% select(ID,Breed) %>% distinct() %>% group_by(Breed) %>% count() %>% filter(n>1) %>% pull(Breed)

# add CBC results
raw <- raw %>% left_join(cbc %>% select(Sample_ID,WBC,RBC))

# Annotate first sample with CBC
best <- raw %>% filter(!is.na(WBC)) 
keep <- best %>% filter(First_Sample)
# for IDs without a First_Sample, use the earliest sample with CBC
keep <- best %>% group_by(ID) %>% summarize(Date_of_sample=min(Date_of_sample)) %>% inner_join(best) %>% bind_rows(keep)
raw <- raw %>% mutate(First_Sample_w_CBC=if_else(Sample_ID %in% keep$Sample_ID,TRUE,FALSE))

# Annotate first sample with weight and age
best <- raw %>% filter(!is.na(Weight_kg)&!is.na(Age_at_sample))
keep <- best %>% filter(First_Sample)
best <- best %>% filter(!(ID %in% keep$ID))
if (length(best$ID)>0){
  keep <- best %>% group_by(ID) %>% summarize(Date_of_sample=min(Date_of_sample)) %>% inner_join(best) %>% bind_rows(keep)
}

raw <- raw %>% mutate(First_Sample_w_AgeWeight=if_else(Sample_ID %in% keep$Sample_ID,TRUE,FALSE))

raw <- raw %>% filter(is.na(replicate)|replicate %in% c("Rep1","High"))
raw <- raw %>% ungroup() %>% mutate(log_DNA_conc=log10(DNA_conc)) %>% select(-DNA_conc) 

raw <- raw %>% filter(!is.na(Date_of_sample)) %>% group_by(ID) %>% summarize(Date_of_first_sample=min(Date_of_sample)) %>% right_join(raw)
raw <- raw %>% mutate(Days_since_first_sample=as.numeric(difftime(Date_of_sample,Date_of_first_sample,units = "days")))

raw <- raw %>% select(ID,Sample_ID,First_Sample,Cohort,Disease_status,Days_since_first_sample,First_Sample_w_AgeWeight,First_Sample_w_CBC,WBC,RBC,Weight_kg,Age_at_sample,log_DNA_conc,TF,Fragment_size_ratio,Plasma_volume_estimated)
raw <- raw %>% pivot_longer(c(log_DNA_conc,TF,Fragment_size_ratio)) %>% rename(stat=name)  %>% filter(!is.na(value)) %>% distinct() 
raw <- raw %>% left_join(dogs)
raw <- raw %>% mutate(Cancer_type3=if_else(Cancer_type=="Lymphoma"|Cancer_type=="Healthy",Cancer_type,"Other cancers"))


# ==============================================================================
# STATISTICAL ANALYSIS (Spearman + Bootstrap)
# ==============================================================================

# Prepare data
pdAll <- raw %>% filter(First_Sample_w_CBC) %>% select(ID,Sample_ID,Cancer_type,Plasma_volume_estimated,WBC,stat,value) %>% left_join(rename_stats)
pdAll <- pdAll %>% mutate(Cancer_type3=if_else(Cancer_type=="Lymphoma","Lymphoma","Other cancers"))
pdAll <- pdAll %>% mutate(Cancer_type3="All cancers") %>% bind_rows(pdAll)
pdAll$Cancer_type3 <- factor(pdAll$Cancer_type3,levels=c("Lymphoma","Other cancers","All cancers"))

# Boot function
boot_spearman <- function(data, indices) {
  d <- data[indices,]
  return(cor(d$value, d$WBC, method = "spearman", use = "complete.obs"))
}

# Calculate Stats
wbc_stats <- pdAll %>% 
  group_by(statistic, stat, Cancer_type3) %>% 
  group_modify(~ {
    d_sub <- .x %>% filter(is.finite(value) & is.finite(WBC))
    
    if(nrow(d_sub) < 3) return(tibble(n=nrow(d_sub), rho=NA, p=NA, ci_low=NA, ci_high=NA))
    
    # Standard Test
    est <- cor.test(d_sub$value, d_sub$WBC, method="spearman", exact=FALSE)
    
    # Bootstrap CI
    set.seed(42)
    boot_res <- boot(data = d_sub, statistic = boot_spearman, R = 2000)
    ci <- tryCatch(boot.ci(boot_res, type="perc")$percent[4:5], error=function(e) c(NA,NA))
    
    tibble(
      n = nrow(d_sub),
      rho = est$estimate,
      p = est$p.value,
      ci_low = ci[1],
      ci_high = ci[2]
    )
  }) %>%
  mutate(
    # Label for Plot
    plot_label = paste0("r = ", round(rho, 2), ", p = ", format_p_plot(p)),
    # Label for CSV
    csv_label = paste0("Spearman: R=", round(rho, 2), " [", round(ci_low, 2), "-", round(ci_high, 2), "]")
  )

write_lines(c("", "Spearman Correlation (WBC vs Metric) with Bootstrapped 95% CI"), outputfile, append=TRUE)
write.table(wbc_stats %>% select(statistic, Cancer_type3, n, rho, p, ci_low, ci_high) %>% mutate(across(c(p, rho, ci_low, ci_high), format_p_csv)), 
            file=outputfile, append=TRUE, sep=",", row.names=FALSE)


# ==============================================================================
# PLOTTING
# ==============================================================================

pdAll <- pdAll %>% mutate(Cancer_color_shape=paste(Cancer_type3,Plasma_volume_estimated))

maxWBC <- ceiling(max(pdAll$WBC, na.rm=T))
colors2 <- c("#ef3b2c","#8073ac")

theme_wbc_corr <- function(){ 
  theme_cowplot(12) %+replace%
    theme(     plot.title = element_text(hjust=0,size=8,face="bold",margin=margin(0,0,4,0)),
               plot.subtitle = element_text(hjust=0,size=7,margin=margin(0,0,4,0),lineheight = 1),
               axis.title = element_text(size=7),axis.text = element_text(size=6.5),
               strip.text.x = element_text(size=6,hjust=0.5,face="bold",color="#f7f7f7"),
               strip.text.y = element_text(size=6,hjust=0.5,face="bold",color="#f7f7f7",angle=90),
               panel.spacing=unit(1.2, "lines"),
               panel.border = element_rect(color = "#525252", fill = NA, size = 1),
               strip.background=element_rect(fill="#525252"),
               legend.position="none",
               legend.justification = c("right", "bottom"),
               legend.text = element_text(size=4.5,vjust=0.5),
               legend.key.size = unit(0.25, 'cm'),
               legend.title = element_blank(),
               legend.margin = margin(-4,1,1,1))
}

## set up dataframe with point shapes and colors
colors <- pdAll %>% select(Cancer_type3) %>% distinct()
colors <- crossing(colors,tibble(Plasma_volume_estimated=c(TRUE,FALSE)))
colors <- colors %>% mutate(color=if_else(Cancer_type3=="Lymphoma","#ef3b2c",if_else(Cancer_type3=="Other cancers","#8073ac","#4d4d4d")))
colors <- colors %>% mutate(shape=if_else(Cancer_type3=="Lymphoma",if_else(Plasma_volume_estimated,21,16),
                                          if_else(Cancer_type3=="Other cancers",if_else(Plasma_volume_estimated,24,17),
                                                  if_else(Plasma_volume_estimated,22,15))))
colors <- colors %>% mutate(Cancer_color_shape=paste(Cancer_type3,Plasma_volume_estimated))
colors <- colors %>% arrange(Cancer_type3,Plasma_volume_estimated) %>% mutate(order=row_number())

# --- Plot Helper Function ---
make_wbc_plot <- function(metric_name, plot_data, stats_data, x_label_scale = NULL) {
  # Filter data
  pd <- plot_data %>% filter(stat == metric_name & !is.na(value))
  pd$Cancer_type3 <- factor(pd$Cancer_type3, levels=c("Lymphoma","Other cancers","All cancers"))
  
  # Setup colors/shapes
  pd$Cancer_color_shape <- factor(pd$Cancer_color_shape, levels=colors %>% arrange(order) %>% pull(Cancer_color_shape))
  shapes <- pd %>% select(Cancer_color_shape) %>% distinct() %>% inner_join(colors) %>% arrange(order) %>% pull(shape)
  
  # Get Stats Labels
  lbls <- stats_data %>% filter(stat == metric_name & Cancer_type3 != "All cancers")
  
  # Assign positions for labels
  # Lymphoma top (maxWBC), Other cancers below it (maxWBC - 3)
  lbls <- lbls %>% 
    mutate(y_pos = if_else(Cancer_type3 == "Lymphoma", maxWBC, maxWBC - 6),
           x_pos = min(pd$value, na.rm=T)) # Left aligned
  
  # Adjust X pos for frag size (log scale)
  if(metric_name == "Fragment_size_ratio") {
    # Since plotting on x=log2(value), labels need log2 x
    lbls$x_pos <- log2(0.25) # Hardcoded to left of plot
  }
  
  # Base Plot
  p <- ggplot(pd, aes(x = if(metric_name=="Fragment_size_ratio") log2(value) else value, 
                      y = WBC, group = Cancer_type3)) + 
    geom_point(aes(color = Cancer_type3, shape = Cancer_color_shape), 
               alpha = 0.5, size = 1, 
               data = pd %>% filter(Cancer_type3 != "All cancers")) +
    geom_smooth(method = 'lm', aes(color = Cancer_type3), 
                linewidth = 0.5, fill = "lightgray", 
                data = pd %>% filter(Cancer_type3 != "All cancers")) +
    
    # Add Text Labels (replacing stat_cor)
    geom_text(data = lbls, aes(x = x_pos, y = y_pos, label = plot_label, color = Cancer_type3),
              hjust = 0, vjust = 1, size = 2.5, show.legend = FALSE) +
    
    scale_color_manual(values = colors2) +
    scale_fill_manual(values = colors2) +
    scale_shape_manual(values = shapes) +
    scale_y_continuous("White blood cell count", limits = c(0, maxWBC)) +
    theme_wbc_corr()
  
  # Axis customization
  if (!is.null(x_label_scale)) {
    p <- p + x_label_scale
  } else {
    p <- p + scale_x_continuous(unique(pd$statistic))
  }
  
  return(p)
}

# --- Generate Plots ---

# 1. Concentration
p_wbc_conc <- make_wbc_plot(
  "log_DNA_conc", pdAll, wbc_stats, 
  scale_x_continuous("cfDNA concentration (ng/mL)", breaks=c(-1:3), labels=10**c(-1:3))
)

# 2. TF
p_wbc_TF <- make_wbc_plot(
  "TF", pdAll, wbc_stats, 
  scale_x_continuous("Tumor fraction")
)

# 3. Fragment Size
xlabels <- c(0.25,0.5,1,2,4)
p_wbc_frag <- make_wbc_plot(
  "Fragment_size_ratio", pdAll, wbc_stats, 
  scale_x_continuous("Fragment size ratio", breaks=log2(xlabels), labels=xlabels)
)

# --- Assemble & Save ---
pdfname <- paste(pdfdir,"Fig_4.pdf",sep="")
p_row_wbc <- plot_grid(p_wbc_conc, p_wbc_TF, p_wbc_frag, ncol=3, label_size=12, labels=LETTERS)

grid <- plot_grid(p_row_wbc, ncol=1, label_size=12)
ggsave(grid, filename=pdfname, width=6.5, height=2.5)
print(paste("Made", pdfname))