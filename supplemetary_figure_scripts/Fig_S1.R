# =========================
# SDOGS — figures: only p= (Wilcoxon for continuous; Fisher for categorical)
# =========================
library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(lubridate)
# This package automatically finds the root folder of the downloaded project.
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ==============================================================================
# 1. SETUP & FILE PATHS
# ==============================================================================

# Define the input files
infile <- here("data", "BB_Supplementary_Data_1.txt")
input_cbc <- here("data", "SCBC.txt")

# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

outputfile <- paste(pdfdir,"Table.Fig_S1.stats.csv",sep="")
write_lines(c("Results for analyses presented in Fig. S1"), outputfile, append=FALSE)

# -------------------------
# Helpers
# -------------------------
.fmt_p <- function(p, digits = 2) {
  out <- formatC(signif(p, digits), digits = digits, format = "g")
  out[is.na(p)] <- NA_character_
  out
}

# Permutation p for unpaired difference in means (continuous)
.perm_p_indep <- function(x, g, B = 10000, seed = 1) {
  set.seed(seed)
  g <- as.factor(g)
  if (nlevels(g) != 2) stop("Group must have 2 levels.")
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(g[ok])
  d_obs <- diff(tapply(x, g, mean, na.rm = TRUE)) # group2 - group1
  perm_stats <- replicate(B, {
    gp <- sample(g)
    # Handle cases where sample might result in empty group (rare but possible in small N)
    means <- tapply(x, gp, mean, na.rm = TRUE)
    if(length(means) < 2) return(FALSE) 
    d  <- diff(means)
    abs(d) >= abs(d_obs)
  })
  (sum(perm_stats) + 1) / (B + 1)
}

# Wrapper for Wilcoxon and Permutation Test
check_wilcox_independent <- function(x, g, B_perm = 10000, seed = 1) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(as.factor(g[ok]))
  if (nlevels(g) != 2) stop("Group factor must have exactly 2 levels.")
  
  # Wilcoxon Rank Sum Test
  wc <- wilcox.test(x ~ g)
  
  # Permuted T-Test (Difference in Means)
  p_perm <- .perm_p_indep(x, g, B = B_perm, seed = seed)
  
  list(
    wilcox_p = wc$p.value,
    perm_t_p = p_perm
  )
}

# =========================
# Data wrangling
# =========================
raw <- as_tibble(read.csv(
  infile,
  header=TRUE,
  sep="\t",
  na.strings = c("NA","N/A","","UNK","?")
))

raw <- raw %>% filter(Sample_ID != "")
raw$X <- NULL
raw <- raw %>% mutate(Sample_type = str_remove(Sample_type,"_plasma"))

# format Date_of_sample
raw <- raw %>%
  mutate(
    Date_of_sample = str_replace(Date_of_sample," \\(approx\\)",""),
    Date_of_sample = parse_date_time(Date_of_sample,"m/d/y")
  )

# lymphoma flag
raw <- raw %>% mutate(Lymphoma = if_else(Cancer_type=="Lymphoma", TRUE, FALSE))

# ancestry
raw <- raw %>%
  filter(!is.na(Breed)) %>%
  mutate(Ancestry = if_else(Breed=="mixed breed","mixed","single")) %>%
  select(ID, Ancestry) %>% distinct() %>%
  right_join(raw, by="ID")

# cancer status
raw <- raw %>%
  mutate(
    Cancer_status     = if_else(Cancer_type!="Healthy","cancer","healthy"),
    Blood_draw_site = if_else(!is.na(Blood_draw_site),
                              paste(Blood_draw_site,Blood_draw_site_side),
                              Blood_draw_site)
  )

# recode spay/sex
raw <- raw %>%
  mutate(
    Spay.neuter.status = if_else(
      Spay.neuter.status=="Y","Yes",
      if_else(Spay.neuter.status=="N","No","other")
    ),
    Sex = if_else(Sex=="F","Female",
                  if_else(Sex=="M","Male","other"))
  )

# timespan (per ID)
raw <- raw %>%
  filter(!is.na(Date_of_sample)) %>%
  group_by(ID) %>%
  summarize(timespan.days = as.numeric(difftime(max(Date_of_sample),
                                                min(Date_of_sample),
                                                units="days")),
            .groups = "drop") %>%
  right_join(raw, by="ID")

# For dogs with Sample_time_point (and missing timespan.days), fallback using that
tmp <- raw %>%
  filter(is.na(timespan.days) & !is.na(Sample_time_point)) %>%
  select(ID, Sample_time_point) %>%
  distinct()

tmp <- tmp %>%
  mutate(
    day = str_replace(str_remove_all(Sample_time_point,"\\."),"Hour","\\."),
    day = floor(as.numeric(str_remove_all(day,"[a-zA-Z ]")))
  ) %>%
  group_by(ID) %>%
  summarize(
    timespan.days2 = if (all(is.na(day))) NA_real_
    else max(day, na.rm=TRUE) - min(day, na.rm=TRUE),
    .groups = "drop"
  )

raw <- raw %>%
  full_join(tmp, by="ID") %>%
  mutate(
    timespan.days = if_else(is.na(timespan.days) & !is.na(timespan.days2),
                            timespan.days2,
                            timespan.days)
  ) %>%
  select(-timespan.days2)

raw <- raw %>%
  mutate(
    timespan.days = if_else(timespan.days==0, 1, timespan.days),
    timespan.days = if_else(ID=="2-BBX", 1, timespan.days)
  )

cbc <- as_tibble(read.csv(input_cbc,
  header=TRUE,
  sep="\t",
  na.strings = c("NA","N/A","","UNK")
))

cbc <- cbc %>% select(Sample_ID) %>% mutate(CBC=TRUE)
cbc <- raw %>% select(ID,Sample_ID) %>% distinct() %>% full_join(cbc, by="Sample_ID")
cbc <- cbc %>% replace_na(list(CBC=FALSE))
raw <- raw %>% left_join(cbc %>% select(Sample_ID, CBC), by="Sample_ID")

d <- raw

# categorize by sampling scheme
tmp <- d %>% filter(timespan.days>=10) %>% pull(timespan.days)
longstr <- paste("longitudinal: ",min(tmp),"-",max(tmp)," days",sep="")

d <- d %>%
  select(ID,Plasma_ID) %>%
  distinct() %>%
  group_by(ID) %>%
  count() %>%
  rename(ndraws = n) %>%
  full_join(d, by="ID")

collection.sets <- d %>%
  filter(timespan.days>=10) %>%
  select(ID) %>% distinct() %>%
  mutate(set = longstr)

collection.sets <- d %>%
  filter(timespan.days>1 & timespan.days<10) %>%
  select(ID) %>% distinct() %>%
  mutate(set="longitudinal: short") %>%
  bind_rows(collection.sets)

collection.sets <- d %>%
  filter(timespan.days<=1 & !is.na(Blood_draw_site)) %>%
  select(ID, Blood_draw_site) %>% distinct() %>%
  group_by(ID) %>% count() %>%
  filter(n>1) %>% select(-n) %>%
  mutate(set="preanalytical: blood draw site") %>%
  bind_rows(collection.sets)

collection.sets <- d %>%
  filter(timespan.days<=1 & !is.na(Time_of_Day) & !(ID %in% collection.sets$ID)) %>%
  select(ID, Time_of_Day) %>% distinct() %>%
  group_by(ID) %>% count() %>%
  filter(n>1) %>% select(-n) %>%
  mutate(set="preanalytical: time of day") %>%
  bind_rows(collection.sets)

collection.sets <- d %>%
  filter(timespan.days<=1 & !is.na(Sample_type)) %>%
  select(ID, Sample_type) %>% distinct() %>%
  group_by(ID) %>% count() %>%
  filter(n>1) %>% select(-n) %>%
  mutate(set="preanalytical: tube type pairs") %>%
  bind_rows(collection.sets)

collection.sets <- d %>%
  filter(timespan.days<=1) %>%
  filter(ndraws>1 & !ID %in% collection.sets$ID) %>%
  select(ID) %>% distinct() %>%
  mutate(set="preanalytical: other") %>%
  bind_rows(collection.sets)

collection.sets <- d %>%
  filter(timespan.days<=1 & ndraws==1 & !(ID %in% collection.sets$ID)) %>%
  select(ID) %>% distinct() %>%
  mutate(set="single sample") %>%
  bind_rows(collection.sets)

collection.sets <- d %>%
  select(ID) %>% distinct() %>%
  full_join(collection.sets, by="ID")

collection.sets <- d %>%
  filter(!is.na(replicate)) %>%
  select(ID) %>% distinct() %>%
  mutate(set="technical replicates") %>%
  bind_rows(collection.sets)

variable_input_dna <- d %>%
  filter(!is.na(ULP_input_DNA))  %>%
  select(ID,Plasma_ID,ULP_input_DNA) %>% distinct() %>%
  group_by(ID,Plasma_ID) %>% count() %>%
  filter(n>1) %>%
  inner_join(d %>% select(ID,Plasma_ID,ULP_input_DNA),
             by=c("ID","Plasma_ID"))

inputDNAstr <- paste(
  "variable input DNA: ",
  min(variable_input_dna$ULP_input_DNA),
  "-",
  max(variable_input_dna$ULP_input_DNA),
  " ng",
  sep=""
)

collection.sets <- d %>%
  filter(ID %in% variable_input_dna$ID) %>%
  select(ID) %>%
  mutate(set=inputDNAstr) %>%
  bind_rows(collection.sets)

collection.sets <- collection.sets %>%
  rename(collection.set = set) %>%
  ungroup()

tube.sets <- d %>%
  filter(Sample_type=="EDTA") %>%
  select(ID) %>% distinct() %>%
  mutate(tube.set="EDTA")

tube.sets <- d %>%
  filter(Sample_type=="Streck") %>%
  select(ID) %>% distinct() %>%
  mutate(tube.set="Streck") %>%
  bind_rows(tube.sets)

cbc.sets <- d %>%
  filter(CBC) %>%
  select(ID,Date_of_sample) %>% distinct() %>%
  group_by(ID) %>% count() %>%
  rename(ncbc=n)

cbc.sets <- d  %>%
  select(ID) %>% distinct() %>%
  full_join(cbc.sets, by="ID") %>%
  replace_na(list(ncbc=0))

cbc.sets <- cbc.sets %>%
  mutate(
    cbc.set = if_else(
      ncbc==0, "no CBC",
      if_else(ncbc==1,"CBC at one timepoint","CBC at more than one timepoint")
    )
  ) %>%
  select(ID,cbc.set) %>%
  distinct() %>%
  ungroup()

# meta about dogs
# Use First_Sample column to identify the canonical metadata row for each dog
dogs <- d %>%
  arrange(ID, desc(First_Sample)) %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup() %>%
  select(ID,Sex,Spay.neuter.status,Ancestry,Cancer_status,Cancer_type,Breed)

dogs <- d %>%
  filter(!is.na(Date_of_sample) & !is.na(TF) & !is.na(DNA_conc)) %>%
  select(ID,Date_of_sample) %>% distinct() %>%
  group_by(ID) %>%
  summarise(Date_first_sample=min(Date_of_sample), .groups="drop") %>%
  right_join(dogs, by="ID")

dogs <- d %>%
  filter(!is.na(Age_at_sample) & !is.na(TF) & !is.na(DNA_conc)) %>%
  select(ID,Age_at_sample) %>% distinct() %>%
  group_by(ID) %>%
  summarise(Age_first_sample=min(Age_at_sample), .groups="drop") %>%
  right_join(dogs, by="ID")

dogs <- raw %>%
  mutate(P30=if_else(str_detect(Cohort,"P30"),TRUE,FALSE)) %>%
  select(ID,P30) %>% distinct() %>%
  # FORCE UNIQUE ID for P30 table before joining
  group_by(ID) %>% slice(1) %>% ungroup() %>%
  right_join(dogs, by="ID")

# first Weight/BCS
weight <- d %>%
  group_by(ID) %>% count() %>%
  filter(n==1) %>%
  inner_join(d, by="ID") %>%
  select(ID,Date_of_sample,Weight_kg) %>%
  distinct()

weight <- d %>%
  filter(!(ID %in% weight$ID) & !is.na(Weight_kg) & !is.na(Date_of_sample)) %>%
  group_by(ID) %>%
  summarize(Date_of_sample = min(Date_of_sample), .groups="drop") %>%
  inner_join(d, by=c("ID","Date_of_sample")) %>%
  select(ID,Date_of_sample,Weight_kg)  %>%
  distinct() %>%
  bind_rows(weight)

# FORCE UNIQUE ID for Weight table before joining
weight <- weight %>% group_by(ID) %>% slice(1) %>% ungroup()

dogs <- dogs %>%
  left_join(weight %>%
              rename(Weight_first_sample = Weight_kg) %>%
              select(-Date_of_sample),
            by="ID") %>%
  distinct()

bcs <- d %>%
  group_by(ID) %>% count() %>%
  filter(n==1) %>%
  inner_join(d, by="ID") %>%
  select(ID,Date_of_sample,BCS) %>%
  distinct()

bcs <- d %>%
  filter(!(ID %in% bcs$ID) & !is.na(BCS) & !is.na(Date_of_sample)) %>%
  group_by(ID) %>%
  summarize(Date_of_sample=min(Date_of_sample), .groups="drop") %>%
  inner_join(d, by=c("ID","Date_of_sample")) %>%
  select(ID,Date_of_sample,BCS)  %>%
  distinct() %>%
  bind_rows(bcs)

# FORCE UNIQUE ID for BCS table before joining
bcs <- bcs %>% group_by(ID) %>% slice(1) %>% ungroup()

dogs <- dogs %>%
  left_join(bcs %>%
              rename(BCS_first_sample = BCS) %>%
              select(-Date_of_sample),
            by="ID") %>%
  distinct()

# FINAL SAFETY CHECK
dogs <- dogs %>% group_by(ID) %>% slice(1) %>% ungroup()

# seq sets
raw <- raw %>%
  mutate(ULPWGS_data = if_else(ULPWGS_data_cf3 | ULPWGS_data_cf4, TRUE, FALSE))

sequencing.sets <- raw %>%
  select(ID,Plasma_ID,ULPWGS_data,Normal_data,Tumor_data,WES_data) %>%
  distinct() %>%
  pivot_longer(c(ULPWGS_data,Normal_data,Tumor_data,WES_data)) %>%
  rename(fnd=value) %>%
  mutate(name = str_remove(name,"_data"))

sequencing.sets <- sequencing.sets %>%
  filter(fnd) %>%
  select(-Plasma_ID) %>%
  distinct() %>%
  arrange(ID,name) %>%
  pivot_wider(names_from=name,values_from=fnd)

sequencing.sets <- sequencing.sets %>%
  full_join(raw %>% select(ID) %>% distinct(), by="ID") %>%
  replace_na(list(ULPWGS=FALSE,WES=FALSE,Normal=FALSE,Tumor=FALSE))

sequencing.sets2 <- sequencing.sets %>%
  filter(ULPWGS) %>%
  mutate(set="Ultra lowpass (ULP)") %>%
  select(ID,set)

sequencing.sets2 <- sequencing.sets %>%
  filter(ULPWGS & Normal & Tumor) %>%
  mutate(set="ULP + paired tumor/normal(WGS)") %>%
  select(ID,set) %>%
  bind_rows(sequencing.sets2)

sequencing.sets2 <- sequencing.sets %>%
  filter(ULPWGS & Normal & WES & Tumor) %>%
  mutate(set="ULP + paired tumor/normal(WGS) + tumor (WES)") %>%
  select(ID,set) %>%
  bind_rows(sequencing.sets2)

sequencing.sets <- sequencing.sets2 %>%
  rename(sequencing.set=set) %>%
  mutate(facet="sequencing\ndata types")

highTF <- raw %>%
  filter(!is.na(TF) & TF>=0.1) %>%
  select(ID) %>% distinct() %>%
  pull(ID)

sequencing.sets <- sequencing.sets %>%
  filter(sequencing.set=="ULP + paired tumor/normal(WGS)") %>%
  filter(ID %in% highTF) %>%
  mutate(sequencing.set=paste(sequencing.set,"(tumor fraction > 10%)")) %>%
  bind_rows(sequencing.sets)

cancer.sets <- d %>%
  select(ID,Cancer_type) %>%
  distinct() %>%
  group_by(Cancer_type) %>% count() %>%
  mutate(cancer.set = if_else(
    Cancer_type=="Healthy" | (n>3 & Cancer_type!="Other"),
    Cancer_type,
    "Other"
  )) %>%
  select(Cancer_type,cancer.set) %>%
  distinct()

cancer.sets <- d %>%
  select(ID,Cancer_type) %>%
  distinct() %>%
  full_join(cancer.sets, by="Cancer_type") %>%
  select(-Cancer_type)

pd <- sequencing.sets %>% rename(set=sequencing.set)
pd <- collection.sets %>% rename(set=collection.set) %>% mutate(facet="collection sets") %>% bind_rows(pd)
pd <- tube.sets %>% rename(set=tube.set) %>% mutate(facet="tube\ntypes") %>% bind_rows(pd)
pd <- cbc.sets %>% rename(set=cbc.set) %>% mutate(facet="CBC") %>% bind_rows(pd)

pd <- pd %>%
  left_join(cancer.sets, by="ID") %>%
  select(ID,cancer.set) %>%
  rename(fill.set=cancer.set) %>%
  distinct() %>%
  full_join(pd, by="ID")

pd <- pd %>%
  select(ID,fill.set,set,facet) %>%
  distinct() %>%
  group_by(fill.set,set,facet) %>%
  count() %>%
  rename(ndogs=n)

pdTotals <- pd %>%
  group_by(facet,set) %>%
  summarize(set.total.n = sum(ndogs), .groups="drop")

levels <- pdTotals %>%
  group_by(facet,set) %>%
  summarize(set.total.n = max(set.total.n), .groups="drop")

levels1 <- c(
  "single sample",
  "preanalytical: blood draw site",
  "preanalytical: time of day",
  "preanalytical: tube type pairs",
  "preanalytical: other",
  "longitudinal: short",
  longstr,
  "technical replicates",
  inputDNAstr,
  "EDTA",
  "Streck",
  "no CBC",
  "CBC at one timepoint",
  "CBC at more than one timepoint"
)

levels2 <- c(
  "Ultra lowpass (ULP)",
  "ULP + normal (WGS)",
  "ULP + tumor (WGS)",
  "ULP + tumor (WGS) (tumor fraction > 10%)",
  "ULP + paired tumor/normal(WGS)",
  "ULP + paired tumor/normal(WGS) (tumor fraction > 10%)",
  "ULP + paired tumor/normal(WGS) + tumor (WES)"
)

levels <- tibble(set=c(rev(levels1),rev(levels2))) %>%
  mutate(order=row_number())

fill.levels <- pd %>%
  filter(fill.set != "Healthy" & facet=="collection sets") %>%
  group_by(fill.set) %>%
  summarize(total=sum(ndogs), .groups="drop")

fill.levels <- rev(c(fill.levels %>% arrange(total) %>% pull(fill.set),"Healthy"))

pd$set      <- factor(pd$set,levels=levels$set)
pd$fill.set <- factor(pd$fill.set,levels=fill.levels)

labels <- pd %>%
  ungroup() %>%
  select(facet,set) %>%
  distinct() %>%
  mutate(sumgrp=row_number())

labels <- pd %>% left_join(labels, by=c("facet","set"))

labels <- tibble(fill.set=rev(fill.levels),
                 fill.order=seq_along(fill.levels)) %>%
  right_join(labels, by="fill.set")

labels <- labels %>%
  arrange(sumgrp,fill.order) %>%
  group_by(sumgrp) %>%
  mutate(xpos=cumsum(ndogs))

labels <- labels %>%
  ungroup() %>%
  rename(label=ndogs) %>%
  rename(ndogs=xpos) %>%
  select(set,fill.set,ndogs,label,facet)

labels <- labels %>%
  filter(label>=5) %>%
  mutate(xpos=((ndogs-label)+ndogs)/2)

samples_per_dog <- raw %>%
  select(ID,Plasma_ID) %>% distinct() %>%
  group_by(ID) %>% count() %>%
  full_join(collection.sets %>% select(ID,collection.set), by="ID")

samples_per_dog <- samples_per_dog %>%
  ungroup() %>%
  group_by(collection.set) %>%
  summarize(min=min(n), max=max(n), .groups="drop")

samples_per_dog <- samples_per_dog %>%
  mutate(
    range = if_else(
      min==max,
      if_else(min==1,paste(min,"sample/dog"),paste(min,"samples/dog")),
      paste(min,"-",max," samples/dog",sep="")
    )
  ) %>%
  rename(set=collection.set) %>%
  select(set,range)

ylabels <- pd %>%
  ungroup() %>%
  select(set) %>% distinct() %>%
  left_join(samples_per_dog, by="set")

ylabels <- ylabels %>%
  mutate(
    ylabel = if_else(
      is.na(range) | str_detect(set,"EDTA") | str_detect(set,"Streck"),
      set,
      paste(set," (",range,")",sep="")
    )
  )

# --------------------------
# Sets overview figure
# --------------------------
title <- ""
p <- ggplot(pd,aes(x=ndogs,y=set)) +
  geom_bar(aes(fill=fill.set),width=0.6,stat="identity",position="stack")

p <- p +
  facet_grid(
    fct_relevel(facet,"sequencing\ndata types","collection sets","tube\ntypes","CBC")~.,
    scales="free_y",
    space="free"
  )

p <- p +
  geom_text(aes(label=set.total.n,x=set.total.n+1),
            hjust=0,size=2.5,data=pdTotals) +
  geom_text(aes(label=label,x=xpos),
            hjust=0.5,size=2,color="#FFFFFF",data=labels)

p <- p +
  scale_fill_manual(values=c("#525252","#d95f02","#1b9e77",
                             "#7570b3","#e7298a","#e6ab02",
                             "#a6761d","#66a61e")) +
  scale_x_continuous("# dogs",limits=c(0,155),breaks=c(0:10)*20) +
  scale_y_discrete("",expand=c(0,0.5),
                   breaks=ylabels$set,
                   labels=ylabels$ylabel) +
  ggtitle(title)

p <- p +
  theme_cowplot(12) +
  theme(
    plot.title    = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size=9),
    axis.text.x  = element_text(size=8),
    axis.text.y  = element_text(size=7),
    legend.position = c(0.98,0.98),
    legend.text  = element_text(size=6),
    legend.title = element_blank(),
    legend.key.size = unit(0.25, 'cm'),
    legend.justification.inside = c(1.05, 1.25),
    legend.box.just = "right",
    legend.margin = margin(2, 3, 3, 3),
    legend.background = element_rect(color="#969696"),
    strip.text.y = element_text(size=6,hjust=0.5,face="bold",color="#f7f7f7"),
    strip.background=element_rect(fill="#525252")
  )

p_nsets_samplesets <- p

table <- pd %>%
  rename("Cancer"=fill.set) %>%
  arrange(facet,set,desc(ndogs))

write_lines(c("","A. Dog counts by sample set"),outputfile,append=TRUE)
write.table(table,file=outputfile,append=TRUE,sep=",",row.names=FALSE)

# ======================================
# Categorical diagnostics (Fisher only)
# ======================================
cat_diag <- bind_rows(
  {
    pdx <- dogs %>% filter(!is.na(Sex), !is.na(Cancer_status))
    tab <- table(pdx$Cancer_status, pdx$Sex)
    ft  <- fisher.test(tab, conf.int = TRUE)
    tibble(
      variable    = "Sex",
      n_healthy   = sum(pdx$Cancer_status=="healthy"),
      n_cancer    = sum(pdx$Cancer_status=="cancer"),
      fisher_p    = ft$p.value,
      odds_ratio  = unname(ft$estimate),
      or_ci_low   = ft$conf.int[1],
      or_ci_high  = ft$conf.int[2],
      min_expected= suppressWarnings(min(chisq.test(tab, correct=FALSE)$expected))
    )
  },
  {
    pdx <- dogs %>% filter(!is.na(Spay.neuter.status), !is.na(Cancer_status))
    tab <- table(pdx$Cancer_status, pdx$Spay.neuter.status)
    ft  <- fisher.test(tab, conf.int = TRUE)
    tibble(
      variable    = "Spay.neuter.status",
      n_healthy   = sum(pdx$Cancer_status=="healthy"),
      n_cancer    = sum(pdx$Cancer_status=="cancer"),
      fisher_p    = ft$p.value,
      odds_ratio  = unname(ft$estimate),
      or_ci_low   = ft$conf.int[1],
      or_ci_high  = ft$conf.int[2],
      min_expected= suppressWarnings(min(chisq.test(tab, correct=FALSE)$expected))
    )
  },
  {
    pdx <- dogs %>% filter(!is.na(Ancestry), !is.na(Cancer_status))
    tab <- table(pdx$Cancer_status, pdx$Ancestry)
    ft  <- fisher.test(tab, conf.int = TRUE)
    tibble(
      variable    = "Ancestry",
      n_healthy   = sum(pdx$Cancer_status=="healthy"),
      n_cancer    = sum(pdx$Cancer_status=="cancer"),
      fisher_p    = ft$p.value,
      odds_ratio  = unname(ft$estimate),
      or_ci_low   = ft$conf.int[1],
      or_ci_high  = ft$conf.int[2],
      min_expected= suppressWarnings(min(chisq.test(tab, correct=FALSE)$expected))
    )
  }
)

write_lines(c("", "B. Fisher exact diagnostics (categorical)"),
            outputfile, append = TRUE)
write.table(cat_diag, file = outputfile, append = TRUE,
            sep = ",", row.names = FALSE)

# =========================
# Plots — categorical panels (figure shows only Fisher p as p=)
# =========================

## Sex
title <- "Sex"
pd_sex <- dogs %>%
  filter(!is.na(Sex), !is.na(Cancer_status)) %>%
  select(ID, Sex, Cancer_status)

tab_sex    <- table(pd_sex$Cancer_status, pd_sex$Sex)
sex_fisher <- fisher.test(tab_sex, conf.int = TRUE)
subtitle   <- paste0("p=", .fmt_p(sex_fisher$p.value))

pdx <- pd_sex %>%
  group_by(Cancer_status) %>%
  count() %>%
  mutate(xlabel = paste(Cancer_status,"\n(N=",n,")",sep="")) %>%
  right_join(pd_sex, by="Cancer_status")

pdx <- pdx %>%
  arrange(Sex) %>%
  group_by(Cancer_status,Sex,xlabel) %>%
  summarize(ndogs=n(), .groups="drop")

pdx <- pdx %>%
  arrange(desc(Sex)) %>%
  group_by(Cancer_status) %>%
  mutate(cum=cumsum(ndogs))

pdx_tot <- pdx %>%
  group_by(Cancer_status) %>%
  summarize(total=sum(ndogs), .groups="drop")

pdx <- pdx %>%
  left_join(pdx_tot, by="Cancer_status") %>%
  mutate(
    frac    = ndogs/total,
    ypos    = cum/total,
    nlabel = paste("N=",ndogs,sep="")
  )

p_sex <- ggplot(pdx,aes(x=xlabel,y=ndogs)) +
  geom_bar(aes(fill=Sex),width=0.5,position="fill",stat="identity") +
  geom_text(aes(y=ypos,label=nlabel),
            size=2.5,vjust=1,nudge_y=-0.01) +
  scale_fill_manual(values=c("#6baed6","#fdbf6f")) +
  scale_y_continuous("Fraction of dogs",breaks=c(0,0.5,1)) +
  ggtitle(title,subtitle=subtitle) +
  theme_cowplot(12) +
  theme(
    plot.title    = element_text(hjust=0,size=10,face="bold"),
    legend.position = "top",
    axis.title.x = element_blank(),
    legend.text  = element_text(size=8),
    legend.title = element_blank(),
    axis.title    = element_text(size=9),
    axis.text     = element_text(size=8)
  )

## Spay/Neuter
title <- "Spay/Neuter status"
pd_spay <- dogs %>%
  filter(!is.na(Spay.neuter.status), !is.na(Cancer_status)) %>%
  select(ID, Spay.neuter.status, Cancer_status)

tab_spay    <- table(pd_spay$Cancer_status, pd_spay$Spay.neuter.status)
spay_fisher <- fisher.test(tab_spay, conf.int = TRUE)
subtitle    <- paste0("p=", .fmt_p(spay_fisher$p.value))

pdx <- pd_spay %>%
  group_by(Cancer_status) %>%
  count() %>%
  mutate(xlabel = paste(Cancer_status,"\n(N=",n,")",sep="")) %>%
  right_join(pd_spay, by="Cancer_status")

pdx <- pdx %>%
  arrange(Spay.neuter.status) %>%
  group_by(Cancer_status,Spay.neuter.status,xlabel) %>%
  summarize(ndogs=n(), .groups="drop")

pdx <- pdx %>%
  arrange(desc(Spay.neuter.status)) %>%
  group_by(Cancer_status) %>%
  mutate(cum=cumsum(ndogs))

pdx_tot <- pdx %>%
  group_by(Cancer_status) %>%
  summarize(total=sum(ndogs), .groups="drop")

pdx <- pdx %>%
  left_join(pdx_tot, by="Cancer_status") %>%
  mutate(
    frac    = ndogs/total,
    ypos    = cum/total,
    nlabel = paste("N=",ndogs,sep="")
  )

p_spay <- ggplot(pdx,aes(x=xlabel,y=ndogs)) +
  geom_bar(aes(fill=Spay.neuter.status),
           width=0.5,position="fill",stat="identity") +
  geom_text(aes(y=ypos,label=nlabel),
            size=2.5,vjust=1,nudge_y=-0.01) +
  scale_fill_manual(values=c("#6baed6","#fdbf6f")) +
  scale_y_continuous("Fraction of dogs",breaks=c(0,0.5,1)) +
  ggtitle(title,subtitle=subtitle) +
  theme_cowplot(12) +
  theme(
    plot.title    = element_text(hjust=0,size=10,face="bold"),
    legend.position = "top",
    axis.title.x = element_blank(),
    legend.text  = element_text(size=8),
    legend.title = element_blank(),
    axis.title    = element_text(size=9),
    axis.text     = element_text(size=8)
  )

## Breed ancestry
title <- "Breed ancestry"
pd_breed <- dogs %>%
  filter(!is.na(Ancestry), !is.na(Cancer_status)) %>%
  select(ID, Ancestry, Cancer_status)

tab_breed    <- table(pd_breed$Cancer_status, pd_breed$Ancestry)
breed_fisher <- fisher.test(tab_breed, conf.int = TRUE)
subtitle     <- paste0("p=", .fmt_p(breed_fisher$p.value))

pdx <- pd_breed %>%
  group_by(Cancer_status) %>%
  count() %>%
  mutate(xlabel = paste(Cancer_status,"\n(N=",n,")",sep="")) %>%
  right_join(pd_breed, by="Cancer_status")

pdx <- pdx %>%
  arrange(Ancestry) %>%
  group_by(Cancer_status,Ancestry,xlabel) %>%
  summarize(ndogs=n(), .groups="drop")

pdx <- pdx %>%
  arrange(desc(Ancestry)) %>%
  group_by(Cancer_status) %>%
  mutate(cum=cumsum(ndogs))

pdx_tot <- pdx %>%
  group_by(Cancer_status) %>%
  summarize(total=sum(ndogs), .groups="drop")

pdx <- pdx %>%
  left_join(pdx_tot, by="Cancer_status") %>%
  mutate(
    frac    = ndogs/total,
    ypos    = cum/total,
    nlabel = paste("N=",ndogs,sep="")
  )

p_breed <- ggplot(pdx,aes(x=xlabel,y=ndogs)) +
  geom_bar(aes(fill=Ancestry),
           width=0.5,position="fill",stat="identity") +
  geom_text(aes(y=ypos,label=nlabel),
            size=2.5,vjust=1,nudge_y=-0.01) +
  scale_fill_manual(values=c("#6baed6","#fdbf6f")) +
  scale_y_continuous("Fraction of dogs",breaks=c(0,0.5,1)) +
  ggtitle(title,subtitle=subtitle) +
  theme_cowplot(12) +
  theme(
    plot.title    = element_text(hjust=0,size=10,face="bold"),
    legend.position = "top",
    axis.title.x = element_blank(),
    legend.text  = element_text(size=8),
    legend.title = element_blank(),
    axis.title    = element_text(size=9),
    axis.text     = element_text(size=8)
  )

# ======================================
# Continuous diagnostics -> CSV; Figure shows Wilcoxon p as p=
# ======================================

# Weight Stats
pdw <- dogs %>% filter(!is.na(Weight_first_sample), !is.na(Cancer_status))
w_stats <- check_wilcox_independent(pdw$Weight_first_sample, pdw$Cancer_status, B_perm=10000)

# BCS Stats
pdb <- dogs %>% filter(!is.na(BCS_first_sample), !is.na(Cancer_status))
b_stats <- check_wilcox_independent(pdb$BCS_first_sample, pdb$Cancer_status, B_perm=10000)

# Age Stats
pda <- dogs %>% filter(!is.na(Age_first_sample), !is.na(Cancer_status))
a_stats <- check_wilcox_independent(pda$Age_first_sample, pda$Cancer_status, B_perm=10000)

# FDR Correction
raw_ps <- c(w_stats$wilcox_p, b_stats$wilcox_p, a_stats$wilcox_p)
fdr_ps <- p.adjust(raw_ps, method = "fdr")

# Create Summary Table
cont_diag_tbl <- tibble(
  variable = c("Weight (kg)", "Body condition score", "Age (years)"),
  n_healthy = c(sum(pdw$Cancer_status=="healthy"), sum(pdb$Cancer_status=="healthy"), sum(pda$Cancer_status=="healthy")),
  n_cancer  = c(sum(pdw$Cancer_status=="cancer"), sum(pdb$Cancer_status=="cancer"), sum(pda$Cancer_status=="cancer")),
  wilcox_p = raw_ps,
  perm_t_p = c(w_stats$perm_t_p, b_stats$perm_t_p, a_stats$perm_t_p),
  fdr_wilcox_p = fdr_ps
)

write_lines(c("", "C. Continuous diagnostics (Wilcoxon, Permutation T-test, FDR)"), outputfile, append = TRUE)
write.table(cont_diag_tbl, file = outputfile, append = TRUE, sep = ",", row.names = FALSE)

# PLOTS (Using raw Wilcoxon p for subtitle)

# Weight Plot
pdx <- pdw %>%
  group_by(Cancer_status) %>%
  count() %>%
  mutate(xlabel = paste(Cancer_status,"\n(N=",n,")",sep="")) %>%
  right_join(pdw, by="Cancer_status")

subtitle <- paste0("p=", .fmt_p(w_stats$wilcox_p))

p_weight <- ggplot(pdx,aes(x=xlabel,y=Weight_first_sample)) +
  geom_boxplot() +
  scale_x_discrete("Cancer status") +
  scale_y_continuous("Weight (kg)",limits=c(0,70)) +
  ggtitle("Weight",subtitle=subtitle) +
  theme_cowplot(12) +
  theme(
    plot.title    = element_text(hjust=0,size=10,face="bold"),
    legend.position = "top",
    axis.title.x = element_blank(),
    legend.text  = element_text(size=8),
    legend.title = element_blank(),
    axis.title    = element_text(size=9),
    axis.text     = element_text(size=8)
  )

# BCS Plot
pdx <- pdb %>%
  group_by(Cancer_status) %>%
  count() %>%
  mutate(xlabel = paste(Cancer_status,"\n(N=",n,")",sep="")) %>%
  right_join(pdb, by="Cancer_status")

subtitle <- paste0("p=", .fmt_p(b_stats$wilcox_p))

p_bcs <- ggplot(pdx,aes(x=xlabel,y=BCS_first_sample)) +
  geom_boxplot() +
  scale_x_discrete("Cancer status") +
  scale_y_continuous("Body condition score",limits=c(0,10)) +
  ggtitle("Body condition score",subtitle=subtitle) +
  theme_cowplot(12) +
  theme(
    plot.title    = element_text(hjust=0,size=10,face="bold"),
    legend.position = "top",
    axis.title.x = element_blank(),
    legend.text  = element_text(size=8),
    legend.title = element_blank(),
    axis.title    = element_text(size=9),
    axis.text     = element_text(size=8)
  )

# Age Plot
pdx <- pda %>%
  group_by(Cancer_status) %>%
  count() %>%
  mutate(xlabel = paste(Cancer_status,"\n(N=",n,")",sep="")) %>%
  right_join(pda, by="Cancer_status")

subtitle <- paste0("p=", .fmt_p(a_stats$wilcox_p))

p_age <- ggplot(pdx,aes(x=xlabel,y=Age_first_sample)) +
  geom_boxplot() +
  scale_x_discrete("Cancer status") +
  scale_y_continuous("Age (years)",limits=c(0,17)) +
  ggtitle("Age", subtitle=subtitle) +
  theme_cowplot(12) +
  theme(
    plot.title    = element_text(hjust=0,size=10,face="bold"),
    legend.position = "top",
    axis.title.x = element_blank(),
    legend.text  = element_text(size=8),
    legend.title = element_blank(),
    axis.title    = element_text(size=9),
    axis.text     = element_text(size=8)
  )

# =========================
# Export figure grid (categoricals + continuous)
# =========================
pdfname <- paste(pdfdir,"Fig_S1.pdf",sep="")
grid <- plot_grid(
  p_sex,p_spay,p_breed,
  p_weight,p_bcs,p_age,
  ncol=3,label_size = 12,labels = "AUTO"
)
ggsave(grid,filename=pdfname,width=6.5,height=5)