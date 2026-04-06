library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(car)
library(broom)
library(lme4)      # Added for Mixed Models
library(lmerTest)  # Added for p-values in Mixed Models
library(lubridate) # Added for date parsing
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

outputfile <- paste(pdfdir,"Table.Fig_3.stats.csv",sep="")
if (!dir.exists(pdfdir)) dir.create(pdfdir, recursive = TRUE)

write_lines(c("Statistical results for analyses presented in Fig. 3 (ANOVA)"),outputfile,append=F)
rename_stats <- tibble(order=c(2,1,3),
                       statistic=c("Tumor fraction","cfDNA concentration (ng/mL)","Fragment size ratio"),
                       stat=c("TF","log_DNA_conc","Fragment_size_ratio")) %>%
  arrange(order)

# Initialize anovaOut as an empty tibble to prevent NULL errors later
anovaOut <- tibble(
  setIn = numeric(),
  modelName = character(),
  sampleNumbers = character(),
  stat = character(),
  Effect = character(),
  DFn = numeric(),
  DFd = numeric(),
  F = numeric(),
  p = numeric(),
  ges = numeric(),
  p_perm = numeric(),
  p_perm_trim = numeric()
)

## ------------------------------------------------------------------
## Permutation helper: Supports Mixed Models with Robust Error Handling
## ------------------------------------------------------------------

nperm_anova <- 100      # number of permutations
trim_frac   <- 0.05     # fraction trimmed

perm_for_model <- function(dsubset, res, formula, nperm = 100, trim_frac = 0.05) {
  out_list <- list()
  stats <- unique(dsubset$stat)
  
  # Check if formula implies a mixed model
  is_mixed <- grepl("\\|", as.character(formula)[3])
  
  for (s in stats) {
    df <- dsubset %>% filter(stat == s)
    if (nrow(df) < 5) next
    
    # Effects we actually care about (Skip ID here)
    effs <- res %>% filter(stat == s, Effect != "ID") %>% pull(Effect)
    if (length(effs) == 0) next
    
    # --- Get Observed F ---
    an_tbl <- tryCatch({
      if (is_mixed) {
        m_obs <- suppressMessages(lmer(formula, data = df, control = lmerControl(optimizer = "bobyqa")))
        at <- suppressMessages(anova(m_obs))
        at <- as.data.frame(at) %>% rownames_to_column("Effect") %>% as_tibble()
        at %>% rename(F = `F value`)
      } else {
        anova_test(data = df, formula = formula) %>% as_tibble()
      }
    }, error = function(e) { NULL })
    
    if(is.null(an_tbl)) next
    
    F_obs  <- an_tbl$F[match(effs, an_tbl$Effect)]
    k      <- length(effs)
    ge_full <- numeric(k)
    ge_trim <- numeric(k)
    
    # --- Determine trimmed subset ---
    n <- nrow(df)
    keep_trim <- integer(0)
    if (n >= 5) {
      q_low  <- quantile(df$value, trim_frac, na.rm = TRUE)
      q_high <- quantile(df$value, 1 - trim_frac, na.rm = TRUE)
      keep_trim <- which(df$value >= q_low & df$value <= q_high)
      min_obs <- if(is_mixed) 4 else length(coef(lm(formula, data=df))) + 1
      if (length(keep_trim) <= min_obs) keep_trim <- integer(0)
    }
    
    if (k > 0) {
      for (i in seq_len(nperm)) {
        df_perm <- df
        df_perm$value <- sample(df$value)
        
        # --- Run Permuted Model ---
        F_perm_full <- rep(NA, k)
        
        if (is_mixed) {
          tryCatch({
            m_perm <- suppressMessages(lmer(formula, data = df_perm, control = lmerControl(optimizer = "bobyqa")))
            if(!is.null(m_perm)) {
              an_perm <- suppressMessages(anova(m_perm))
              an_perm <- as.data.frame(an_perm) %>% rownames_to_column("Effect")
              F_perm_full <- an_perm$`F value`[match(effs, an_perm$Effect)]
            }
          }, error = function(e) { NULL })
        } else {
          tryCatch({
            an_perm_full <- anova_test(data = df_perm, formula = formula) %>% as_tibble()
            F_perm_full  <- an_perm_full$F[match(effs, an_perm_full$Effect)]
          }, error = function(e) { NULL })
        }
        
        ge_full <- ge_full + as.numeric(!is.na(F_perm_full) & !is.na(F_obs) & F_perm_full >= F_obs)
        
        # --- Run Trimmed Permuted Model ---
        if (length(keep_trim) > 0) {
          df_perm_trim <- df_perm[keep_trim, , drop = FALSE]
          F_perm_trim <- rep(NA, k)
          
          if (is_mixed) {
            tryCatch({
              m_perm_trim <- suppressMessages(lmer(formula, data = df_perm_trim, control = lmerControl(optimizer = "bobyqa")))
              if(!is.null(m_perm_trim)) {
                an_perm_trim <- suppressMessages(anova(m_perm_trim))
                an_perm_trim <- as.data.frame(an_perm_trim) %>% rownames_to_column("Effect")
                F_perm_trim <- an_perm_trim$`F value`[match(effs, an_perm_trim$Effect)]
              }
            }, error = function(e) { NULL })
          } else {
            tryCatch({
              an_perm_trim <- anova_test(data = df_perm_trim, formula = formula) %>% as_tibble()
              F_perm_trim  <- an_perm_trim$F[match(effs, an_perm_trim$Effect)]
            }, error = function(e) { NULL })
          }
          ge_trim <- ge_trim + as.numeric(!is.na(F_perm_trim) & !is.na(F_obs) & F_perm_trim >= F_obs)
        }
      }
    }
    
    p_perm      <- if (k > 0) ge_full / nperm else numeric(0)
    p_perm_trim <- if (length(keep_trim) > 0 && k > 0) ge_trim / nperm else rep(NA_real_, k)
    
    out_list[[length(out_list) + 1]] <- tibble(
      stat        = s,
      Effect      = effs,
      p_perm      = p_perm,
      p_perm_trim = p_perm_trim
    )
  }
  
  if (length(out_list) == 0) {
    tibble(stat = character(0), Effect = character(0), p_perm = double(0), p_perm_trim = double(0))
  } else {
    bind_rows(out_list)
  }
}

## ------------------------------------------------------------------
## Data Processing
## ------------------------------------------------------------------

raw <- as_tibble(read.csv(infile,header=T,sep="\t",na.strings = c("NA","N/A","","UNK")))
cbc <- as_tibble(read.csv(infile_cbc,header=T,sep="\t",na.strings = c("NA","N/A","","UNK")))

raw <- raw %>% mutate(Sample_type=str_remove(Sample_type,"_plasma"))

raw <- raw %>%
  mutate(Date_of_sample=str_replace(Date_of_sample," \\(approx\\)","")) %>%
  mutate(Date_of_sample=parse_date_time(Date_of_sample,"m/d/y"))

raw <- raw %>%
  mutate(Sex.spay=if_else(is.na(Spay.neuter.status),
                          Spay.neuter.status,
                          paste(Sex,Spay.neuter.status,sep=".")))
dogs <-  raw %>% select(ID,Cancer_type,Sex,Sex.spay,Breed) %>% distinct()
raw <- raw %>% select(-Cancer_type,-Sex,-Breed,-Sex.spay,-Spay.neuter.status) %>% distinct()

dogs <- dogs %>%
  filter(!is.na(Breed)) %>%
  mutate(Ancestry=if_else(Breed=="mixed breed","mixed","single")) %>%
  select(ID,Ancestry) %>%
  distinct() %>%
  right_join(dogs, by = "ID")

breeds <- dogs %>%
  select(ID,Breed) %>%
  distinct() %>%
  group_by(Breed) %>%
  count() %>%
  mutate(Breed_v2=if_else(n>=2,Breed,"other breeds")) %>%
  rename(ndogs_in_breed=n) %>%
  select(Breed,ndogs_in_breed,Breed_v2) %>%
  distinct()
dogs <- dogs %>% left_join(breeds, by = "Breed")

dogs <- dogs %>% mutate(Cancer_status=if_else(Cancer_type!="Healthy","cancer","healthy"))
dogs <- dogs %>% mutate(Lymphoma=if_else(Cancer_type=="Lymphoma",TRUE,FALSE))

okbreeds_gt1 <- dogs %>%
  select(ID,Breed) %>%
  distinct() %>%
  group_by(Breed) %>%
  count() %>%
  filter(n>1) %>%
  pull(Breed)

raw <- raw %>% left_join(cbc %>% select(Sample_ID,WBC,RBC), by = "Sample_ID")

best <- raw %>% filter(!is.na(WBC))
keep <- best %>% filter(First_Sample)
keep <- best %>%
  group_by(ID) %>%
  summarize(Date_of_sample=min(Date_of_sample), .groups = "drop") %>%
  inner_join(best, by = c("ID","Date_of_sample")) %>%
  bind_rows(keep)
raw <- raw %>% mutate(First_Sample_w_CBC=if_else(Sample_ID %in% keep$Sample_ID,TRUE,FALSE))

best <- raw %>% filter(!is.na(Weight_kg)&!is.na(Age_at_sample))
keep <- best %>% filter(First_Sample)
best <- best %>% filter(!(ID %in% keep$ID))
if (length(best$ID)>0){
  keep <- best %>%
    group_by(ID) %>%
    summarize(Date_of_sample=min(Date_of_sample), .groups = "drop") %>%
    inner_join(best, by = c("ID","Date_of_sample")) %>%
    bind_rows(keep)
}
raw <- raw %>% mutate(First_Sample_w_AgeWeight=if_else(Sample_ID %in% keep$Sample_ID,TRUE,FALSE))

raw <- raw %>% filter(is.na(replicate)|replicate=="Rep1")
raw <- raw %>% ungroup() %>% mutate(log_DNA_conc=log10(DNA_conc)) %>% select(-DNA_conc)

raw <- raw %>%
  filter(!is.na(Date_of_sample)) %>%
  group_by(ID) %>%
  summarize(Date_of_first_sample=min(Date_of_sample), .groups = "drop") %>%
  right_join(raw, by = "ID")
raw <- raw %>%
  mutate(Days_since_first_sample=as.numeric(difftime(Date_of_sample,Date_of_first_sample,units = "days")))

raw <- raw %>%
  select(ID,Sample_ID,First_Sample,Cohort,Disease_status,Days_since_first_sample,
         First_Sample_w_AgeWeight,First_Sample_w_CBC,WBC,RBC,Weight_kg,Age_at_sample,
         log_DNA_conc,TF,Fragment_size_ratio,Plasma_volume_estimated,BCS,Sequencing_group,Sample_type)
raw <- raw %>%
  pivot_longer(c(log_DNA_conc,TF,Fragment_size_ratio)) %>%
  rename(stat=name)  %>%
  filter(!is.na(value)) %>%
  distinct()
raw <- raw %>% left_join(dogs, by = "ID")
raw <- raw %>%
  mutate(Cancer_type3=if_else(Cancer_type=="Lymphoma"|Cancer_type=="Healthy",
                              Cancer_type,"Other cancers"))

## VIFs
vifs <- raw %>%
  filter(!is.na(Disease_status)&!is.na(Cancer_type)) %>%
  group_by(stat) %>%
  do(tidy(vif(lm(value~Disease_status+Cancer_type, .))))
vifs <- vifs %>% mutate(model="value~Disease_status+Cancer_type")
vifs2 <- raw %>%
  filter(!is.na(Breed)&!is.na(Cancer_type)&(Breed %in% okbreeds_gt1)) %>%
  group_by(stat) %>%
  do(tidy(vif(lm(value~Breed+Cancer_type, .))))
vifs <- vifs2 %>% mutate(model="value~Breed+Cancer_type") %>% bind_rows(vifs)
vifs3 <- raw %>%
  filter(!is.na(Breed)&!is.na(Weight_kg)&(Breed %in% okbreeds_gt1)) %>%
  group_by(stat) %>%
  do(tidy(vif(lm(value~Breed+Weight_kg, .))))
vifs <- vifs3 %>% mutate(model="value~Breed+Weight") %>% bind_rows(vifs)
vifs4 <- raw %>%
  filter(!is.na(Cancer_type)&!is.na(Sequencing_group)) %>%
  group_by(stat) %>%
  do(tidy(vif(lm(value~Cancer_type+Sequencing_group, .))))
vifs <- vifs4 %>% mutate(model="value~Cancer_type+Sequencing_group") %>% bind_rows(vifs)
vifs5 <- raw %>%
  filter(!is.na(Cancer_type)&!is.na(Weight_kg)) %>%
  group_by(stat) %>%
  do(tidy(vif(lm(value~Cancer_type+Weight_kg, .))))
vifs <- vifs5 %>% mutate(model="value~Cancer_type+Weight") %>% bind_rows(vifs)
vifs6 <- raw %>%
  filter(!is.na(BCS)&!is.na(Weight_kg)) %>%
  group_by(stat) %>%
  do(tidy(vif(lm(value~BCS+Weight_kg, .))))
vifs <- vifs6 %>% mutate(model="value~BCS+Weight") %>% bind_rows(vifs)

vifs$GVIF  <- vifs$x[,1]
vifs$Df    <- vifs$x[,2]
vifs$GVIF2 <- vifs$x[,3]
vifs <- vifs %>% rename('GVIF^(1/(2*Df))'=GVIF2 )
vifs <- rename_stats %>% left_join(vifs, by = "stat")
vifs <- vifs %>% select(model,statistic,GVIF,Df,'GVIF^(1/(2*Df))')

writeLetter <- 1
write_lines(c("",paste(LETTERS[writeLetter],". Anova analysis of factors affecting liquid biopsy results",sep="")),
            outputfile,append=TRUE)
write.table(vifs %>% arrange(model,statistic),
            file=outputfile,append=TRUE,sep=",",row.names=FALSE)
writeLetter <- writeLetter + 1

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

#### Make ANOVA panel
d <- raw %>%
  filter(!is.na(value)) %>%
  filter(Cancer_status=="cancer") %>%
  select(-Sex)
d <- d %>% mutate(value=if_else(stat=="Fragment_size_ratio",log2(value),value))

rename_effects <- tibble(
  Effect=rev(c("ID","Cancer_type","Disease_status","Cancer_status","Lymphoma","Osteosarcoma",
               "Breed","WBC","RBC","Days_since_first_sample","Sex.spay","Ancestry",
               "Weight_kg","Age_at_sample","Plasma_volume_estimated","BCS","Sequencing_group", "Sample_type")),
  label=rev(c("Dog ID","Cancer type","Disease status","Healthy or has cancer?","Lymphoma?",
              "Osteosarcoma?","Breed","White blood cell count","Red blood cell count",
              "Days since first sample","Sex","Mixed vs. single breed","Weight","Age",
              "Volume estimated?","Body Condition Score","Sequencing group", "Sample type"))
)
rename_effects <- rename_effects %>% mutate(order=row_number())
nsamplesID <- NULL
nsamples <- NULL

# ==============================================================================
# MODEL 1: Disease_status + (1|ID) [UPDATED TO MIXED MODEL]
# ==============================================================================
modelName <- "Disease_status (Mixed Model)"
dsubset <- d %>%
  filter(Cancer_status=="cancer") %>%
  filter(!is.na(Disease_status))
ntotal <- dsubset %>% group_by(stat) %>% count()
nsamples <- length(unique(dsubset$Sample_ID))
ndogs <- length(unique(dsubset$ID))
sampleNumbers <- paste(ndogs," dogs",sep="")

# Use lmer() for mixed model analysis
res_list <- list()
icc_list <- list()
stats_m1 <- unique(dsubset$stat)

for(s in stats_m1) {
  df_s <- dsubset %>% filter(stat == s)
  
  m_res <- tryCatch({
    m <- lmer(value ~ Disease_status + (1|ID), data = df_s, control = lmerControl(optimizer = "bobyqa"))
    
    # Extract ANOVA table
    an_tbl <- anova(m) %>% as.data.frame() %>% rownames_to_column("Effect") %>% as_tibble()
    an_tbl <- an_tbl %>% 
      rename(F = `F value`, p = `Pr(>F)`) %>%
      mutate(
        stat = s, 
        DFn = NumDF, 
        DFd = DenDF,
        ges = (F * DFn) / (F * DFn + DFd)
      )
    
    # --- CALCULATE ICC FOR DOG ID ---
    vars <- as.data.frame(VarCorr(m))
    var_id <- vars[vars$grp=="ID", "vcov"]
    var_resid <- vars[vars$grp=="Residual", "vcov"]
    icc <- var_id / (var_id + var_resid)
    
    row_id <- tibble(
      Effect = "ID",
      DFn = NA, DFd = NA, F = NA, p = NA, 
      stat = s, 
      ges = icc
    )
    
    list(an_tbl = an_tbl, row_id = row_id)
  }, error = function(e) { print(paste("Model 1 Failed for", s, ":", e$message)); NULL })
  
  if(!is.null(m_res)) {
    res_list[[s]] <- m_res$an_tbl
    icc_list[[s]] <- m_res$row_id
  }
}
res <- bind_rows(res_list)

# Permutation for Mixed Model
perm <- perm_for_model(dsubset, res, value ~ Disease_status + (1|ID),
                       nperm = nperm_anova, trim_frac = trim_frac)
res <- res %>% left_join(perm, by = c("stat","Effect"))

# ADDED BACK: Bind ICC rows so they appear in CSV
res <- bind_rows(res, bind_rows(icc_list))

anovaOut <- res %>%
  mutate(setIn=1,modelName=modelName,sampleNumbers=sampleNumbers) %>%
  left_join(ntotal, by = "stat") %>%
  bind_rows(anovaOut)

# ==============================================================================
# MODEL 3: Demographics (Fixed Effect, First Sample)
# ==============================================================================
modelName <- "Cancer_type+Sex+Age_at_Sample+Breed+Plasma_est+Seq_group"
dsubset <- d %>%
  filter(First_Sample) %>%
  filter(Cancer_status == "cancer") %>%
  filter(
    !is.na(Cancer_type) &
      !is.na(Sex.spay) &
      !is.na(Age_at_sample) &
      !is.na(Breed) &
      !is.na(Plasma_volume_estimated) &
      !is.na(Sequencing_group) &
      !is.na(Sample_type)
  )
dsubset <- dsubset %>% filter(nsamples == 1 | First_Sample)
ntotal <- dsubset %>% group_by(stat) %>% count()
nsamples <- length(unique(dsubset$Sample_ID))
ndogs    <- length(unique(dsubset$ID))
sampleNumbers <- paste(ndogs, " dogs", sep = "")

## 3a) Concentration
d_conc <- dsubset %>% filter(stat == "log_DNA_conc")
preds_conc <- c("Cancer_type", "Sex.spay", "Age_at_sample", "Breed", "Plasma_volume_estimated")
if(n_distinct(d_conc$Sequencing_group, na.rm=TRUE) > 1) preds_conc <- c(preds_conc, "Sequencing_group")
if(n_distinct(d_conc$Sample_type, na.rm=TRUE) > 1) preds_conc <- c(preds_conc, "Sample_type")
form_conc <- as.formula(paste("value ~", paste(preds_conc, collapse = " + ")))

res_conc <- d_conc %>%
  group_by(stat) %>%
  anova_test(form_conc) %>%
  as_tibble()
perm_conc <- perm_for_model(d_conc, res_conc, form_conc, nperm = nperm_anova, trim_frac = trim_frac)
res_conc <- res_conc %>% left_join(perm_conc, by = c("stat","Effect"))

## 3b) Other Stats
d_other <- dsubset %>% filter(stat != "log_DNA_conc")
stats_other <- unique(d_other$stat)
res_other_list <- list()

for(s in stats_other) {
  d_s <- d_other %>% filter(stat == s)
  preds_s <- c("Cancer_type", "Sex.spay", "Age_at_sample", "Breed")
  if(n_distinct(d_s$Sequencing_group, na.rm=TRUE) > 1) preds_s <- c(preds_s, "Sequencing_group")
  if(n_distinct(d_s$Sample_type, na.rm=TRUE) > 1) preds_s <- c(preds_s, "Sample_type")
  form_s <- as.formula(paste("value ~", paste(preds_s, collapse = " + ")))
  
  an_s <- anova_test(data = d_s, formula = form_s) %>% as_tibble()
  an_s <- an_s %>% mutate(stat = s)
  perm_s <- perm_for_model(d_s, an_s, form_s, nperm = nperm_anova, trim_frac = trim_frac)
  an_s <- an_s %>% left_join(perm_s, by = c("stat","Effect"))
  res_other_list[[s]] <- an_s
}
res_other <- bind_rows(res_other_list)
anova_model <- bind_rows(res_conc, res_other)

anovaOut <- anova_model %>%
  mutate(setIn = if_else(nrow(anovaOut) == 0, 1, max(anovaOut$setIn) + 1), modelName = modelName, sampleNumbers = sampleNumbers) %>%
  left_join(ntotal, by = "stat") %>%
  bind_rows(anovaOut)

# ==============================================================================
# MODEL 4: Physical (Fixed Effect, First Sample)
# ==============================================================================
modelName <- "Cancer_type+Weight_kg+BCS+covariates"
dsubset <- d %>%
  filter(First_Sample) %>%
  filter(Cancer_status=="cancer") %>%
  filter(
    !is.na(Cancer_type) &
      !is.na(Weight_kg) &
      !is.na(BCS) & 
      !is.na(Sequencing_group) & 
      !is.na(Sample_type)
  )
dsubset <- dsubset %>% filter(nsamples==1|First_Sample)
ntotal <- dsubset %>% group_by(stat) %>% count()
nsamples <- length(unique(dsubset$Sample_ID))
ndogs <- length(unique(dsubset$ID))
sampleNumbers <- paste(ndogs," dogs",sep="")

stats_4 <- unique(dsubset$stat)
res_4_list <- list()

for(s in stats_4) {
  d_s <- dsubset %>% filter(stat == s)
  preds_s <- c("Cancer_type", "Weight_kg", "BCS")
  if(s == "log_DNA_conc" && n_distinct(d_s$Plasma_volume_estimated, na.rm=TRUE) > 1) {
    preds_s <- c(preds_s, "Plasma_volume_estimated")
  }
  if(n_distinct(d_s$Sequencing_group, na.rm=TRUE) > 1) preds_s <- c(preds_s, "Sequencing_group")
  if(n_distinct(d_s$Sample_type, na.rm=TRUE) > 1) preds_s <- c(preds_s, "Sample_type")
  form_s <- as.formula(paste("value ~", paste(preds_s, collapse = " + ")))
  
  an_s <- anova_test(data = d_s, formula = form_s) %>% as_tibble()
  an_s <- an_s %>% mutate(stat = s)
  perm_s <- perm_for_model(d_s, an_s, form_s, nperm = nperm_anova, trim_frac = trim_frac)
  an_s <- an_s %>% left_join(perm_s, by = c("stat","Effect"))
  res_4_list[[s]] <- an_s
}
res <- bind_rows(res_4_list)

anovaOut <- res %>%
  mutate(setIn=max(anovaOut$setIn)+1,modelName=modelName,sampleNumbers=sampleNumbers) %>%
  left_join(ntotal, by = "stat") %>%
  bind_rows(anovaOut)

# ==============================================================================
# MODEL 5: WBC (Mixed Model) [UPDATED]
# ==============================================================================
modelName <- "WBC + (1|ID)"
dsubset <- d  %>% filter(!is.na(WBC)) %>% filter(Cancer_status=="cancer")
ntotal <- dsubset  %>% group_by(stat) %>% count()
ndogs <- dsubset %>% select(ID) %>% distinct() %>% count()
nsamples <- length(unique(dsubset$Sample_ID))
sampleNumbers <- paste(ndogs," dogs",sep="")

# Use lmer()
res_list <- list()
icc_list <- list()
stats_m5 <- unique(dsubset$stat)

for(s in stats_m5) {
  df_s <- dsubset %>% filter(stat == s)
  m_res <- tryCatch({
    m <- lmer(value ~ WBC + (1|ID), data = df_s, control = lmerControl(optimizer = "bobyqa"))
    an_tbl <- anova(m) %>% as.data.frame() %>% rownames_to_column("Effect") %>% as_tibble()
    an_tbl <- an_tbl %>% 
      rename(F = `F value`, p = `Pr(>F)`) %>%
      mutate(
        stat = s, 
        DFn = NumDF, 
        DFd = DenDF, 
        ges = (F * DFn) / (F * DFn + DFd)
      )
    
    vars <- as.data.frame(VarCorr(m))
    var_id <- vars[vars$grp=="ID", "vcov"]
    var_resid <- vars[vars$grp=="Residual", "vcov"]
    icc <- var_id / (var_id + var_resid)
    
    row_id <- tibble(
      Effect = "ID",
      DFn = NA, DFd = NA, F = NA, p = NA, 
      stat = s, 
      ges = icc
    )
    list(an_tbl = an_tbl, row_id = row_id)
  }, error = function(e) { print(paste("Model 5 Failed for", s, ":", e$message)); NULL })
  
  if(!is.null(m_res)) {
    res_list[[s]] <- m_res$an_tbl
    icc_list[[s]] <- m_res$row_id
  }
}
res <- bind_rows(res_list)

perm <- perm_for_model(dsubset, res, value ~ WBC + (1|ID),
                       nperm = nperm_anova, trim_frac = trim_frac)
res <- res %>% left_join(perm, by = c("stat","Effect"))

# Bind ICC rows so they appear in CSV
res <- bind_rows(res, bind_rows(icc_list))

anovaOut <- res %>%
  mutate(setIn=max(anovaOut$setIn)+1,modelName=modelName,sampleNumbers=sampleNumbers) %>%
  left_join(ntotal, by = "stat") %>%
  bind_rows(anovaOut)

# ==============================================================================
# MODEL 6: Cancer type + WBC (Fixed Effect, First Sample with CBC)
# ==============================================================================
modelName <- "Cancer type+WBC"
dsubset <- d  %>% filter(First_Sample_w_CBC) %>% filter(Cancer_status=="cancer")
ntotal <- dsubset %>% group_by(stat) %>% count()
ndogs <- dsubset %>% select(ID) %>% distinct() %>% count()
nsamples <- length(unique(dsubset$Sample_ID))
sampleNumbers <- paste(ndogs," dogs",sep="")

res <- dsubset %>%
  group_by(stat) %>%
  anova_test(value ~ Cancer_type + WBC) %>%
  as_tibble()
perm <- perm_for_model(dsubset, res, value ~ Cancer_type + WBC,
                       nperm = nperm_anova, trim_frac = trim_frac)
res <- res %>% left_join(perm, by = c("stat","Effect"))

anovaOut <- res %>%
  mutate(setIn=max(anovaOut$setIn)+1,modelName=modelName,sampleNumbers=sampleNumbers) %>%
  left_join(ntotal, by = "stat") %>%
  bind_rows(anovaOut)

### Write Consolidated ANOVA results to output file
# Combine everything into one table with all p-values
final_stats <- anovaOut %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "fdr")) %>%
  mutate(sig = case_when(
    is.na(p) ~ "not significant",
    p <= 0.05 ~ "significant",
    TRUE ~ "not significant"
  )) %>%
  select(
    setIn, modelName, stat, Effect,
    DFn, DFd, `F`, p, p.adj,
    p_perm, p_perm_trim,
    ges, sig
  ) %>%
  arrange(setIn, stat, desc(ges))

write_lines(c("", "Combined Statistical Results (Raw, FDR, Permuted, Robust)"), outputfile, append=TRUE)
write.table(final_stats, file=outputfile, append=TRUE, sep=",", row.names=FALSE)

anovaOut <- anovaOut %>% mutate(set=paste("#",setIn,") ",sampleNumbers,sep=""))

### Facet labels etc.
anovaOut <- anovaOut %>% left_join(rename_stats %>% select(-order), by = "stat")
anovaOut$statistic <- factor(anovaOut$statistic,
                             levels=(rename_stats %>% arrange(order) %>% pull(statistic)))

facetrow_levels <- anovaOut %>%
  select(set,setIn) %>%
  distinct() %>%
  arrange(setIn) %>%
  pull(set)
anovaOut$set <- factor(anovaOut$set,levels=facetrow_levels)

rename_effects <- rename_effects %>% arrange(order)
# Use unique() to prevent duplicating "ID" which is already in rename_effects
anovaOut$Effect <- factor(anovaOut$Effect, levels=unique(c(rename_effects$Effect)))

anovaOut <- anovaOut %>%
  mutate(pstr=if_else(p<0.001,
                      paste("p=",format(p,digits=1,scientific=TRUE),sep=""),
                      paste("p=",round(p,3),sep="")))

### Make ANOVA plot
# EXPLICITLY REMOVE "ID" (ICC) FROM THE PLOT DATA HERE
pd <- anovaOut %>%
  mutate(sig = case_when(
    is.na(p) ~ "not significant",
    p <= 0.05 ~ "significant",
    TRUE ~ "not significant"
  )) %>%
  filter(Effect != "ID") %>% 
  select(setIn,modelName,Effect,ges,sig,pstr,set,statistic,stat,n)

sample_counts <- pd %>%
  select(set,statistic,n,Effect) %>%
  distinct() %>%
  mutate(nstr=paste("N=",n,sep=""))

# Ensure N is shown by attaching it to the primary effect of each model
sample_counts <- sample_counts %>%
  filter(Effect %in% c("Disease_status", "Cancer_type", "WBC"))

rename_effects <- rename_effects %>%
  mutate(label=if_else(label=="Weight","Weight*",label))
rename_effects <- rename_effects %>%
  mutate(label=if_else(label=="Body Condition Score","Body Condition Score*",label))

p <- ggplot(pd,aes(x=ges,y=Effect))
p <- p + geom_bar(aes(fill=sig),stat="identity",width=0.5)
p <- p + geom_vline(xintercept = 0,color="#969696",linewidth=0.2)
p <- p + geom_text(aes(label=pstr,color=sig),hjust=0,size=1.75,nudge_x = 0.01)
p <- p + geom_text(aes(label=nstr),x=1.12,hjust=1,vjust=0.5,
                   color="#852119",data=sample_counts,size=1.75)
p <- p + facet_grid(rows=vars(set),cols=vars(statistic),space="free_y",scales="free_y")
# Manually add ID to labels
p <- p + scale_y_discrete("", breaks=rename_effects$Effect, labels=rename_effects$label)
p <- p + scale_color_manual(values=c("#525252","#000000"))
p <- p + scale_fill_manual(values=c("#969696","#000000"))
p <- p + scale_x_continuous("Anova effect size (ges)",
                            limits=c(0,1.1),
                            expand = expansion(mult = c(0, .05)),
                            breaks=c(0,0.25,0.50,0.75,1))
p <- p + theme_anova()
p_anova <- p

#### Boxplots and PDF assembly (unchanged)
theme_compare <- function(){
  theme_cowplot(12) %+replace%
    theme(     plot.title = element_text(hjust=0,size=8,face="bold",margin=margin(0,0,4,0)),
               plot.subtitle = element_text(hjust=0,size=7,margin=margin(0,0,4,0),lineheight = 1),
               axis.title.x=element_blank(),
               axis.title.y = element_text(size=8,angle=90,vjust=1),axis.text.x = element_text(size=5),
               axis.text.y = element_text(size=6),
               axis.line=element_line(linewidth=0.25),axis.ticks=element_line(linewidth=0.25),
               legend.position="none")
}

pdAll <- raw %>% ungroup() %>% filter(First_Sample) %>% mutate(Cancer_type4=Cancer_type3)
pdAll <- pdAll %>% filter(Cancer_type!="Healthy") %>%
  mutate(Cancer_type4="All cancers") %>%
  bind_rows(pdAll)
pdAll <- pdAll %>% select(ID,Sample_ID,Plasma_volume_estimated,Cancer_type3,Cancer_type4,stat,value)
pdAll <- pdAll %>% mutate(Cancer_status=if_else(Cancer_type4=="Healthy","Healthy","Cancer"))
pdAll <- pdAll %>% left_join(rename_stats, by = "stat")
pdAll$Cancer_type4 <- factor(pdAll$Cancer_type4,
                             levels=c("Healthy","All cancers","Lymphoma","Other cancers"))
pdAll <- pdAll %>%
  mutate(Cancer_color_shape=paste(Cancer_type3,Plasma_volume_estimated))

# Colors and Shapes logic for Filled (Estimated) vs Open (Measured)
colors <- pdAll %>% select(Cancer_type3) %>% distinct()
colors <- crossing(colors,tibble(Plasma_volume_estimated=c(TRUE,FALSE)))
colors <- colors %>%
  mutate(color=if_else(Cancer_type3=="Lymphoma","#ef3b2c",
                       if_else(Cancer_type3=="Other cancers","#8073ac","#4d4d4d")))

# Explicitly map shape based on User request:
# TRUE (Estimated) -> Filled (16/17/15)
# FALSE (Measured) -> Open (1/2/0)
colors <- colors %>%
  mutate(shape = case_when(
    Cancer_type3 == "Lymphoma" & Plasma_volume_estimated ~ 16, # Filled Circle
    Cancer_type3 == "Lymphoma" & !Plasma_volume_estimated ~ 1, # Open Circle
    Cancer_type3 == "Other cancers" & Plasma_volume_estimated ~ 17, # Filled Triangle
    Cancer_type3 == "Other cancers" & !Plasma_volume_estimated ~ 2, # Open Triangle
    TRUE & Plasma_volume_estimated ~ 15, # Filled Square (Healthy/All)
    TRUE & !Plasma_volume_estimated ~ 0  # Open Square
  )) %>%
  mutate(Cancer_color_shape=paste(Cancer_type3,Plasma_volume_estimated)) %>%
  arrange(Cancer_type3,Plasma_volume_estimated) %>%
  mutate(order=row_number())

# Create named vectors for robust scale mapping
shape_map <- setNames(colors$shape, colors$Cancer_color_shape)
color_map <- setNames(colors$color, colors$Cancer_type3)

# Shape Map for Panel B & C (Cancer Type Only)
# Healthy=Square (15), Lymphoma=Circle (16), Other=Triangle (17)
shape_map_simple <- c("Healthy"=15, "Lymphoma"=16, "Other cancers"=17)

# Define Base Colors for manual scale
colors3 <- c("Healthy"="#525252", "Lymphoma"="#ef3b2c", "Other cancers"="#8073ac")

pd <- pdAll %>% filter(stat=="log_DNA_conc") %>% filter(!is.na(value))
pd$Cancer_color_shape <- factor(pd$Cancer_color_shape, levels=colors$Cancer_color_shape)

xlabels <- pd %>%
  group_by(Cancer_type4) %>%
  count() %>%
  mutate(xlabel=paste(Cancer_type4,"\n(N=",n,")",sep="")) %>%
  mutate(xlabel=str_replace(xlabel,"\\_","\n"))

p <- ggplot(pd,aes(x=Cancer_type4,y=value))
p <- p + geom_boxplot(color="#525252",outlier.size=1,outlier.shape=NA)
p <- p + geom_jitter(aes(color=Cancer_type3, shape=Cancer_color_shape),
                     size=1, alpha=0.5, height=0, width=0.2)
p <- p + stat_compare_means(method = "t.test",
                            comparisons= list(
                              c("Healthy","All cancers"),
                              c("Lymphoma","Other cancers"),
                              c("Healthy","Lymphoma"),
                              c("Healthy","Other cancers")),
                            size=2.5)
p <- p + scale_color_manual(values=colors3)
p <- p + scale_shape_manual(values=shape_map)
p <- p + scale_x_discrete("Cancer type",
                          breaks=xlabels$Cancer_type4,
                          labels=xlabels$xlabel)
p <- p + scale_y_continuous(unique(pd$statistic),
                            breaks=c(-1:3),
                            labels=10**c(-1:3),
                            expand = expansion(mult = c(0.05, .1)))
p <- p + theme_compare()
p_type_conc <- p

pd <- pdAll %>% filter(stat=="TF") %>% filter(!is.na(value))
# Use Simple Shape for TF (Panel B)
pd$Cancer_color_shape <- factor(pd$Cancer_color_shape, levels=colors$Cancer_color_shape)

xlabels <- pd %>%
  group_by(Cancer_type4) %>%
  count() %>%
  mutate(xlabel=paste(Cancer_type4,"\n(N=",n,")",sep="")) %>%
  mutate(xlabel=str_replace(xlabel,"\\_","\n"))

p <- ggplot(pd,aes(x=Cancer_type4,y=value))
p <- p + geom_boxplot(color="#525252",outlier.size=1,outlier.shape=NA)
# Use Cancer_type3 for shape directly
p <- p + geom_jitter(aes(color=Cancer_type3, shape=Cancer_type3),
                     size=1, alpha=0.5, height=0, width=0.2)
p <- p + stat_compare_means(method = "t.test",
                            comparisons= list(
                              c("Healthy","All cancers"),
                              c("Lymphoma","Other cancers"),
                              c("Healthy","Lymphoma"),
                              c("Healthy","Other cancers")),
                            size=2.5)
p <- p + scale_color_manual(values=colors3)
p <- p + scale_shape_manual(values=shape_map_simple)
p <- p + scale_x_discrete("Cancer type",
                          breaks=xlabels$Cancer_type4,
                          labels=xlabels$xlabel)
p <- p + scale_y_continuous(unique(pd$statistic),
                            expand = expansion(mult = c(0.05, .1)))
p <- p + theme_compare()
p_type_TF <- p

ylabels <- c(0.25,0.5,1,2,4)
pd <- pdAll %>% filter(stat=="Fragment_size_ratio") %>% filter(!is.na(value))
# Use Simple Shape for Frag (Panel C)

xlabels <- pd %>%
  group_by(Cancer_type4) %>%
  count() %>%
  mutate(xlabel=paste(Cancer_type4,"\n(N=",n,")",sep="")) %>%
  mutate(xlabel=str_replace(xlabel,"\\_","\n"))

p <- ggplot(pd,aes(x=Cancer_type4,y=log2(value)))
p <- p + geom_boxplot(color="#525252",outlier.size=1,outlier.shape=NA)
# UPDATED: Use Cancer_type3 for shape directly
p <- p + geom_jitter(aes(color=Cancer_type3, shape=Cancer_type3),
                     size=1, alpha=0.5, height=0, width=0.2)
p <- p + stat_compare_means(method = "t.test",
                            comparisons= list(
                              c("Healthy","All cancers"),
                              c("Lymphoma","Other cancers"),
                              c("Healthy","Lymphoma"),
                              c("Healthy","Other cancers")),
                            size=2.5)
p <- p + scale_color_manual(values=colors3)
p <- p + scale_shape_manual(values=shape_map_simple)
p <- p + scale_x_discrete("Cancer type",
                          breaks=xlabels$Cancer_type4,
                          labels=xlabels$xlabel)
p <- p + scale_y_continuous(unique(pd$statistic),
                            expand = expansion(mult = c(0.05, .1)),
                            breaks=log2(ylabels),labels=ylabels)
p <- p + theme_compare()
p_type_frag <- p

# --- ADDED: LEGEND EXTRACTION & INSERTION ---
# 1. Color Legend (Cancer Type)
p_legend_color_source <- ggplot(pd, aes(x=Cancer_type4, y=value, color=Cancer_type3)) +
  geom_point(size=2.5) + # Smaller size
  scale_color_manual(
    values=colors3, 
    name="Cancer Type", 
    labels=c("Healthy", "Lymphoma", "Other cancers"),
    # Override shapes to show Square/Circle/Triangle
    guide = guide_legend(override.aes = list(shape = c(15, 16, 17), size=2.5))
  ) +
  theme_cowplot(10) + theme(legend.position="bottom", legend.justification="center", legend.margin=margin(0,0,0,0))
legend_c <- get_legend(p_legend_color_source)

# 2. Shape Legend (Plasma Volume)
# Dummy data to force shapes 16 (Filled) vs 1 (Open)
legend_data <- tibble(
  label = factor(c("Estimated", "Measured"), levels = c("Estimated", "Measured")),
  x = c(1, 2), y = c(1, 1)
)

p_legend_shape_source <- ggplot(legend_data, aes(x=x, y=y, shape=label)) +
  geom_point(size=2.5, color="black") +
  scale_shape_manual(name="Plasma Vol.", values=c("Estimated"=16, "Measured"=1)) +
  theme_cowplot(10) + theme(legend.position="bottom", legend.justification="center", legend.margin=margin(0,0,0,0))
legend_s <- get_legend(p_legend_shape_source)

# Combine legends
legend_row <- plot_grid(legend_c, legend_s, nrow=1, rel_widths=c(1, 0.6))

# Assemble plots
pdfname <- paste(pdfdir,"Fig_3.pdf",sep="")
row1 <- plot_grid(p_type_conc, p_type_TF, p_type_frag, labels=LETTERS[1:3], nrow=1, label_size = 12)

# Insert legend between row1 and row2
grid <- plot_grid(row1, legend_row, p_anova, ncol=1, labels=c("", "", LETTERS[4]), label_size=12, rel_heights=c(1, 0.1, 1))

ggsave(grid, filename=pdfname, width=6.5, height=7.5)
print(paste("Made",pdfname))