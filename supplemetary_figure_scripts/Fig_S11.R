library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(car)
library(broom)
library(stats)
library(coin)
library(lubridate)
library(lme4)      # Added for Mixed Models
library(lmerTest)  # Added for p-values in LMM
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
infile_cbc <- here("data", "SCBC.txt")
  
# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

outputfile <- paste(pdfdir,"Table.S11_stats.csv",sep="")
write_lines(c("Statistical results for analyses presented in Fig. S11"),outputfile,append=FALSE)

rename_stats <- tibble(
  order=c(2,1,3),
  statistic=c("Tumor fraction","cfDNA concentration (ng/mL)","Fragment size ratio"),
  stat=c("TF","log_DNA_conc","Fragment_size_ratio") 
) %>% arrange(order)

# --- Helper: Formatter ---
format_p_csv <- function(x) {
  sapply(x, function(val) {
    if (is.na(val) || is.infinite(val)) return("NA")
    if (val == 0) return("< 2.2e-16") 
    if (val < 0.001) {
      formatC(val, format = "e", digits = 4)
    } else {
      formatC(val, format = "f", digits = 4)
    }
  })
}

### Process input files
if(file.exists(infile) & file.exists(infile_cbc)) {
  raw <- as_tibble(read.csv(infile,header=T,sep="\t",na.strings = c("NA","N/A","","UNK")))
  raw <- raw %>% filter(Sample_ID!="")
  if("X" %in% colnames(raw)) raw$X <- NULL
  
  cbc <- as_tibble(read.csv(infile_cbc,header=T,sep="\t",na.strings = c("NA","N/A","","UNK")))
  cbc <- cbc %>% filter(Sample_ID!="")
  
  raw <- raw %>% mutate(Sample_type=str_remove(Sample_type,"_plasma")) %>% filter(Cohort %in% c("LSA_U01_treatment_cohort_1","LSA_U01_treatment_cohort_2"))
  
  # format Date_of_sample column as a date
  raw <- raw %>% mutate(Date_of_sample=str_replace(Date_of_sample," \\(approx\\)","")) %>% mutate(Date_of_sample=parse_date_time(Date_of_sample,"m/d/y"))
  raw <- raw %>% mutate(Sex.spay=if_else(is.na(Spay.neuter.status),Spay.neuter.status,paste(Sex,Spay.neuter.status,sep=".")))
  dogs <-  raw %>% select(ID,Sex,Sex.spay,Breed) %>% distinct()
  raw <- raw %>% select(-Sex,-Breed,-Sex.spay,-Spay.neuter.status) %>% distinct()
  # add column for mixed or single breed
  dogs <- dogs %>% filter(!is.na(Breed)) %>% mutate(Ancestry=if_else(Breed=="mixed breed","mixed","single")) %>% select(ID,Ancestry) %>% distinct() %>% right_join(dogs)
  
  breeds <- dogs %>% select(ID,Breed) %>% distinct() %>% group_by(Breed) %>% count() %>% mutate(Breed_v2=if_else(n>=2,Breed,"other breeds")) %>% rename(ndogs_in_breed=n) %>% select(Breed,ndogs_in_breed,Breed_v2) %>% distinct()
  dogs <- dogs %>% left_join(breeds)
  
  # add CBC results
  raw <- raw %>% left_join(cbc %>% select(Sample_ID,WBC,RBC))
  
  # Annotate first sample with CBC
  best <- raw %>% filter(!is.na(WBC)) 
  keep <- best %>% filter(First_Sample)
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
  
  raw <- raw %>% filter(is.na(replicate)|replicate=="Rep1")
  raw <- raw %>% ungroup() %>% mutate(log_DNA_conc=log10(DNA_conc)) %>% select(-DNA_conc) 
  
  raw <- raw %>% filter(!is.na(Date_of_sample)) %>% group_by(ID) %>% summarize(Date_of_first_sample=min(Date_of_sample)) %>% right_join(raw)
  raw <- raw %>% mutate(Days_since_first_sample=as.numeric(difftime(Date_of_sample,Date_of_first_sample,units = "days")))
  
  raw <- raw %>% select(ID,Sample_ID,Stage,First_Sample,Cohort,Disease_status,Days_since_first_sample,First_Sample_w_AgeWeight,First_Sample_w_CBC,WBC,RBC,Weight_kg,Age_at_sample,log_DNA_conc,TF,Fragment_size_ratio,Plasma_volume_estimated,BCS,Sequencing_group,Tumor_size_cm)
  raw <- raw %>% pivot_longer(c(log_DNA_conc,TF,Fragment_size_ratio)) %>% rename(stat=name)  %>% filter(!is.na(value)) %>% distinct() 
  raw <- raw %>% left_join(dogs)
  
  #### Make ANOVA panel 
  d <- raw %>% filter(!is.na(value)) %>% select(-Sex) 
  # Calculate Ratio
  d$size_weight <- d$Tumor_size_cm/d$Weight_kg
  
  d <- d %>% mutate(value = if_else(stat == "Fragment_size_ratio" & value > 0, log2(value), value)) %>% filter(is.finite(value))
} else {
  print("Warning: Input files not found. Creating dummy 'd' object for structure validation.")
  d <- tibble(
    ID = factor(1:20),
    stat = rep(c("TF", "log_DNA_conc", "Fragment_size_ratio"), length.out=20),
    value = rnorm(20),
    Disease_status = factor(rep(c("A","B"),10)),
    Tumor_size_cm = runif(20, 1, 10),
    Weight_kg = runif(20, 5, 30),
    Age_at_sample = runif(20, 2, 12),
    BCS = factor(rep(1:5, 4)),
    size_weight = runif(20, 0.1, 0.5),
    Stage = factor(rep(1:4, 5)),
    Breed = factor(rep(c("A","B"), 10)),
    First_Sample = TRUE
  )
}

theme_anova <- function(){ 
  theme_cowplot(12) %+replace%
    theme(     plot.title = element_text(hjust=0,size=7,margin=margin(0,0,4,0)),
               plot.subtitle = element_text(hjust=0,size=7,margin=margin(0,0,4,0),lineheight = 1),
               legend.position = "none",
               axis.line.x=element_line(linewidth=0.25),axis.ticks.x=element_line(linewidth=0.25),
               axis.line.y=element_blank(),
               axis.ticks.y=element_blank(),
               legend.text = element_text(size=7),
               legend.title = element_blank(),
               axis.title = element_text(size=6), 
               axis.text.x = element_text(size=5.5),
               axis.text.y = element_text(size=5.5,hjust=1,vjust=0.5),
               strip.text.x = element_text(size=5,hjust=0.5,face="bold",color="#f7f7f7",lineheight=0.9),
               strip.text.y = element_text(size=5,hjust=0,vjust=1,color="#f7f7f7",angle=0),
               panel.spacing=unit(.2, "lines"),
               panel.border = element_rect(color = "#525252", fill = NA, linewidth = 1),
               strip.background=element_rect(fill="#525252"),
               panel.grid.major.x = element_line(color="#bdbdbd",linewidth = 0.15))
}

rename_effects <- tibble(Effect=rev(c("ID","Disease_status","Stage","Cancer_status","Tumor_size_cm","Lymphoma","Osteosarcoma","Breed","WBC","RBC","Days_since_first_sample","Sex.spay","Ancestry","Weight_kg","Age_at_sample","Plasma_volume_estimated","BCS","Sequencing_group","size_weight")),label=rev(c("Dog ID","Disease status","Disease stage","Healthy or has cancer?","Tumor size (cm)","Lymphoma?","Osteosarcoma?","Breed","White blood cell count","Red blood cell count","Days since first sample","Sex","Mixed vs. single breed","Weight","Age","Volume estimated?","Body Condition Score","Sequencing group","Ratio (Size/Weight)")))
rename_effects <- rename_effects %>% mutate(order=row_number())


# ==============================================================================
# HELPER FUNCTION: RUN ANOVA SUITE (LMM or Standard)
# ==============================================================================
run_anova_suite <- function(data, formula_str, set_id, model_name) {
  
  # Detect if this is a Mixed Model
  is_mixed <- grepl("\\|", formula_str)
  f <- as.formula(formula_str)
  
  # Prepare Robust Formula: Remove "ID" terms to prevent saturation in robust check
  clean_form <- formula_str
  if (grepl("ID", formula_str)) {
    clean_form <- gsub("\\+ \\(1\\|ID\\)", "", clean_form) 
    clean_form <- gsub("\\+ ID|ID \\+|ID", "", clean_form)
    clean_form <- gsub("~ \\+", "~", clean_form)
  }
  f_robust <- as.formula(clean_form)
  
  res <- data %>%
    group_by(stat) %>%
    group_modify(~ {
      
      # 1. Main Model (LMM or ANOVA)
      if(is_mixed) {
        # Linear Mixed Model
        model_fit <- tryCatch({
          suppressMessages(lmer(f, data = .x))
        }, error = function(e) return(NULL))
        
        if (is.null(model_fit)) return(tibble())
        
        # Anova Table
        std_anova <- anova(model_fit) %>% 
          as.data.frame() %>% 
          rownames_to_column("Effect") %>% 
          as_tibble() %>%
          rename(F = `F value`, p = `Pr(>F)`) %>%
          mutate(
            DFn = NumDF, 
            DFd = DenDF,
            ges = (F * DFn) / (F * DFn + DFd)
          )
        
      } else {
        # Standard ANOVA
        std_anova <- tryCatch({
          rstatix::anova_test(.x, f) %>% as_tibble()
        }, error = function(e) return(NULL))
        
        if (is.null(std_anova)) return(tibble())
      }
      
      std_anova$p_perm <- NA
      std_anova$p_robust <- NA
      
      # 2. Try Full Robust Model (Fixed Effects Only approximation)
      rob_mod <- NULL
      n_obs <- nrow(.x)
      n_preds <- length(attr(terms(f_robust), "term.labels"))
      
      # Only run FULL robust model if we have enough data to avoid crash
      if(n_preds > 0 && n_obs > (5 * n_preds)) {
        # tryCatch specifically for the "dimnames" error in lmrob
        rob_mod <- tryCatch(
          lmrob(f_robust, data = .x, setting="KS2014"), 
          error=function(e) NULL
        )
      }
      
      # Iterate effects
      for (i in 1:nrow(std_anova)) {
        eff <- std_anova$Effect[i]
        
        # --- ROBUST P-VALUE LOGIC ---
        p_rob_val <- NA
        
        # Strategy A: Extract from Full Model
        if (!is.null(rob_mod)) {
          # tryCatch for summary/coefficient extraction failures
          try({
            coefs <- summary(rob_mod)$coefficients
            matches <- rownames(coefs)[grep(paste0("^", eff), rownames(coefs))]
            if (length(matches) > 0) {
              p_vals <- coefs[matches, "Pr(>|t|)", drop=FALSE]
              p_vals <- p_vals[!is.na(p_vals)]
              if(length(p_vals) > 0) p_rob_val <- min(p_vals)
            }
          }, silent=TRUE)
        }
        
        # Strategy B: Fallback to Univariate Robust Test
        if (is.na(p_rob_val)) {
          f_uni <- as.formula(paste("value ~", eff))
          try({
            # Nested tryCatch for the Univariate model creation AND summary
            p_rob_val <- tryCatch({
              uni_mod <- lmrob(f_uni, data = .x, setting="KS2014")
              coefs <- summary(uni_mod)$coefficients
              if(nrow(coefs) > 1) {
                min(coefs[2:nrow(coefs), "Pr(>|t|)"], na.rm=TRUE)
              } else {
                NA
              }
            }, error = function(e) NA)
          }, silent=TRUE)
        }
        std_anova$p_robust[i] <- p_rob_val
        
        # --- PERMUTATION P-VALUE ---
        f_uni <- as.formula(paste("value ~", eff))
        d_uni <- .x %>% select(value, all_of(eff)) %>% drop_na() %>% as.data.frame() # Force DF
        
        # Skip if constant variance or singular
        run_perm <- TRUE
        if(nrow(d_uni) < 5) run_perm <- FALSE
        if(run_perm) {
          # Check diversity of predictor
          if(length(unique(d_uni[[eff]])) < 2) run_perm <- FALSE
          # Check for constant response (Variance check)
          if(var(d_uni$value) <= 1e-12) run_perm <- FALSE
        }
        
        if(run_perm) {
          if(is.character(d_uni[[eff]]) || (is.numeric(d_uni[[eff]]) && length(unique(d_uni[[eff]])) < 5)) {
            d_uni[[eff]] <- as.factor(d_uni[[eff]])
          }
          
          p_u <- tryCatch({
            # Explicit coin calls
            it <- coin::independence_test(f_uni, data=d_uni, distribution=coin::approximate(nresample=2000))
            coin::pvalue(it)
          }, error=function(e) {
            # message(paste("Permutation fail for", eff, ":", e$message)) # Debug if needed
            NA 
          })
          
          std_anova$p_perm[i] <- as.numeric(p_u)
        }
      }
      return(std_anova)
    }) %>%
    ungroup()
  
  # Add Metadata
  count_data <- data %>% group_by(stat) %>% count()
  ndogs <- length(unique(data$ID))
  sample_str <- paste(ndogs, " dogs", sep="")
  
  res <- res %>% 
    left_join(count_data, by="stat") %>%
    mutate(
      setIn = set_id, 
      modelName = model_name, 
      sampleNumbers = sample_str
    )
  
  return(res)
}


# ==============================================================================
# RUN MODELS
# ==============================================================================

anovaOut <- tibble()

# --- Model 1: Status + Size + Weight + Age + BCS + (1|ID) ---
dsubset <- d %>% 
  filter(!is.na(Disease_status) & !is.na(Tumor_size_cm) & 
           !is.na(Weight_kg) & !is.na(Age_at_sample) & !is.na(BCS))

mod1 <- run_anova_suite(dsubset, 
                        "value ~ Disease_status + Tumor_size_cm + Weight_kg + Age_at_sample + BCS + (1|ID)", 
                        1, "Status + Size + Wt + Age + BCS + (1|ID)")
anovaOut <- bind_rows(anovaOut, mod1)


# --- Model 2: Size + Weight + Ratio + (1|ID) ---
dsubset <- d %>% 
  filter(!is.na(Tumor_size_cm) & !is.na(Weight_kg) & !is.na(size_weight))

mod2 <- run_anova_suite(dsubset, 
                        "value ~ Tumor_size_cm + Weight_kg + size_weight + (1|ID)", 
                        2, "Size + Wt + Ratio + (1|ID)")
anovaOut <- bind_rows(anovaOut, mod2)


# --- Model 3: Stage + Size + Breed + Wt + Age + BCS (FIRST SAMPLES) ---
dsubset <- d %>% 
  filter(First_Sample == TRUE | First_Sample == "TRUE") %>% 
  filter(!is.na(Stage) & !is.na(Tumor_size_cm) & !is.na(Breed) &
           !is.na(Weight_kg) & !is.na(Age_at_sample) & !is.na(BCS))

mod3 <- run_anova_suite(dsubset, 
                        "value ~ Stage + Tumor_size_cm + Breed + Weight_kg + Age_at_sample + BCS", 
                        3, "First: Stage + Size + Breed + Wt + Age + BCS")
anovaOut <- bind_rows(anovaOut, mod3)


# ==============================================================================
# FORMAT OUTPUT
# ==============================================================================

anovaOut <- anovaOut %>% mutate(set=paste("#",setIn,") ",sampleNumbers,sep=""))

# Add FDR Adjustment (Global)
anovaOut <- anovaOut %>% 
  mutate(
    p.adj = p.adjust(p, method="fdr"),
    p.perm.adj = p.adjust(p_perm, method="fdr"),
    p.robust.adj = p.adjust(p_robust, method="fdr")
  )

# Plot Significance based on raw p-value
anovaOut <- anovaOut %>% mutate(sig = if_else(p <= 0.05, "significant", "not significant")) 

# Join Labels
anovaOut <- anovaOut %>% left_join(rename_stats %>% select(-order)) 
anovaOut$statistic <- factor(anovaOut$statistic, levels=(rename_stats %>% arrange(order) %>% pull(statistic)))

facetrow_levels <- anovaOut %>% select(set,setIn) %>% distinct() %>% arrange(setIn) %>% pull(set)
anovaOut$set <- factor(anovaOut$set, levels=facetrow_levels)

rename_effects <- rename_effects %>% arrange(order) 
anovaOut$Effect <- factor(anovaOut$Effect, levels=rename_effects$Effect)

# Plot Labels based on raw p-value
anovaOut <- anovaOut %>% 
  mutate(pstr = if_else(p < 0.001, 
                        paste("p=", format(p, digits=1, scientific=TRUE), sep=""), 
                        paste("p=", round(p, 3), sep="")))

# SAVE CSV (Includes FDR and raw)
write_out <- anovaOut %>% 
  select(set, modelName, statistic, Effect, DFn, DFd, F, ges, p, p.adj, p_perm, p.perm.adj, p_robust, p.robust.adj) %>%
  mutate(across(starts_with("p"), format_p_csv))

write.table(write_out, file=outputfile, append=TRUE, sep=",", row.names=FALSE)


# ==============================================================================
# PLOTTING
# ==============================================================================

pd <- anovaOut
pd <- anovaOut %>% select(setIn, modelName, Effect, ges, sig, pstr, set, statistic, stat, n)

# Counts annotations
sample_counts <- pd %>% select(set, statistic, n, Effect) %>% distinct() %>% mutate(nstr = paste("N=", n, sep=""))
# Attach N counts to primary variable of interest for each model to avoid clutter
sample_counts <- sample_counts %>% 
  filter(Effect %in% c("Disease_status", "Tumor_size_cm", "Stage")) %>%
  group_by(set, statistic) %>%
  slice(1) # Keep one per panel

p <- ggplot(pd, aes(x=ges, y=Effect)) 
p <- p + geom_bar(aes(fill=sig), stat="identity", width=0.5)
p <- p + geom_vline(xintercept = 0, color="#969696", linewidth=0.2)
p <- p + geom_text(aes(label=pstr, color=sig), hjust=0, size=1.75, nudge_x = 0.01)
p <- p + geom_text(aes(label=nstr), x=1.12, hjust=1, vjust=0.5, color="#852119", data=sample_counts, size=1.75)
p <- p + facet_grid(rows=vars(set), cols=vars(statistic), space="free_y", scales="free_y")
p <- p + scale_y_discrete("", breaks=rename_effects$Effect, labels=rename_effects$label)
p <- p + scale_color_manual(values=c("#525252","#000000"))
p <- p + scale_fill_manual(values=c("#969696","#000000"))
p <- p + scale_x_continuous("Anova effect size (ges)", limits=c(0,1.2), expand = expansion(mult = c(0, .05)), breaks=c(0,0.25,0.50,0.75,1), labels=c("0","0.25","0.50","0.75","1"))
p <- p + theme_anova()
p_anova <- p 

pdfname <- paste(pdfdir,"Fig_S11.pdf",sep="")
grid <- plot_grid(p_anova, ncol=1, label_size = 12, labels="", rel_heights = c(1,1))
ggsave(grid, filename=pdfname, width=6.5, height=3.75)
print(paste("Made", pdfname))