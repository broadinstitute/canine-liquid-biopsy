library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(ggrastr)

# This package automatically finds the root folder of the downloaded project.
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ==============================================================================
# 1. SETUP & FILE PATHS
# ==============================================================================

# Define the input file (pointing to the 'data' folder in your repo)
infile <- here("data", "BB_Supplementary_Data_1.txt")
infile_cbc <- here("data", "SCBC.txt")
infile_conmuts <- here("data", "CONMUTS.txt")
infile_depth <- here("data","READDEPTH.txt")
infile_metrics <- here("data", "ULPMETRICS.txt")

# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

pdfname <- paste(pdfdir,"Fig_2.pdf",sep="")

### Make first panel 
raw <- as_tibble(read.csv(infile,header=T,sep="\t",na.strings = c("NA","N/A","","UNK","?")))

metrics <- as_tibble(read.csv(infile_metrics,header=T,sep="\t",na.strings = c("NA","N/A","","UNK","?")))
metrics <- metrics %>% rename(ULP_input_DNA=Total.ng.cfDNA.Input) %>% filter(!is.na(ULP_input_DNA))
metrics <- metrics %>% select(Sample_ID,ULP_input_DNA)  %>% distinct()
raw <- raw %>% left_join(metrics)

# get samples with tumor fraction > 10% (for including in panel B)
samples_over_10per <- raw %>% filter(TF>0.1) %>% select(ID) %>% distinct()

dognames <- raw %>% select(ID,Sample_ID,Plasma_ID) %>% distinct()

raw <- raw %>% mutate(Sample_type=str_remove(Sample_type,"_plasma"))

# format Date_of_sample column as a date
raw <- raw %>% mutate(Date_of_sample=str_replace(Date_of_sample," \\(approx\\)","")) %>% mutate(Date_of_sample=parse_date_time(Date_of_sample,"m/d/y"))

# add column for whether cancer type is lymphoma or something else 
raw <- raw %>% mutate(Lymphoma=if_else(Cancer_type=="Lymphoma",TRUE,FALSE))
raw %>% filter(!is.na(Time_of_Day)&is.na(Blood_draw_site)) 

# add column for mixed or single breed
raw <- raw %>% filter(!is.na(Breed)) %>% mutate(Ancestry=if_else(Breed=="mixed breed","mixed","single")) %>% select(ID,Ancestry) %>% distinct() %>% right_join(raw)

# recode Cancer_status as healthy/cancer
raw <- raw %>% mutate(Cancer_status=if_else(Cancer_type!="Healthy","cancer","healthy"))

raw <- raw %>% mutate(Blood_draw_site=if_else(!is.na(Blood_draw_site),paste(Blood_draw_site,Blood_draw_site_side),Blood_draw_site))

# recode Spay/neuter as yes/no, sex as Female/Male
raw <- raw %>% mutate(Spay.neuter.status=if_else(Spay.neuter.status=="Y","Yes",if_else(Spay.neuter.status=="N","No","other")))
raw <- raw %>% mutate(Sex=if_else(Sex=="F","Female",if_else(Sex=="M","Male","other")))

# get total number of days of sampling
raw <- raw %>% filter(!is.na(Date_of_sample)) %>% group_by(ID) %>% summarize(timespan.days=as.numeric(difftime(max(Date_of_sample),min(Date_of_sample),unit="days"))) %>% right_join(raw)

# get total number of days of sampling for samples without sampling date 
tmp <- raw %>% filter(is.na(timespan.days)&!is.na(Sample_time_point)) %>% select(ID,Sample_time_point) %>% distinct()
tmp <- tmp %>% mutate(day=str_replace(str_remove_all(Sample_time_point,"\\."),"Hour","\\.")) %>% mutate(day=floor(as.numeric(str_remove_all(day,"[a-zA-Z ]"))))
tmp <- tmp %>% group_by(ID) %>% summarize(timespan.days2=max(day)-min(day))
raw <- raw %>% full_join(tmp) %>% mutate(timespan.days=if_else(is.na(timespan.days)&!is.na(timespan.days2),timespan.days2,timespan.days)) %>% select(-timespan.days2)
raw <- raw %>% mutate(timespan.days=if_else(timespan.days==0,1,timespan.days))
raw <- raw %>% mutate(timespan.days=if_else(ID=="2-BBX",1,timespan.days))

cbc <- as_tibble(read.csv(infile_cbc,header=T,sep="\t",na.strings = c("NA","N/A","","UNK")))
cbc <- cbc %>% select(Sample_ID) %>% mutate(CBC=TRUE)
cbc <- raw %>% select(ID,Sample_ID) %>% distinct() %>% full_join(cbc)
cbc <- cbc %>% replace_na(list(CBC=FALSE)) 

raw <- raw %>% left_join(cbc)

d <- raw 

### categorize dogs by type of sampling scheme 

tmp <- d %>% filter(timespan.days>=10) %>% pull(timespan.days)
longstr <- paste("longitudinal: ",min(tmp),"-",max(tmp)," days",sep="")
d <- d %>% select(ID,Plasma_ID) %>% distinct() %>% group_by(ID) %>% count() %>% rename(ndraws=n) %>% full_join(d)
collection.sets <- d %>% filter(timespan.days>=10) %>% select(ID) %>% distinct() %>% mutate(set=longstr)
collection.sets <- d %>% filter(timespan.days>1&timespan.days<10) %>% select(ID) %>% distinct() %>% mutate(set="longitudinal: short") %>% bind_rows(collection.sets)
collection.sets <- d %>% filter(timespan.days<=1&!is.na(Blood_draw_site)) %>% select(ID,Blood_draw_site) %>% distinct() %>% group_by(ID) %>% count() %>% filter(n>1) %>% select(-n) %>% mutate(set="preanalytical: blood draw site") %>% bind_rows(collection.sets)
collection.sets <- d %>% filter(timespan.days<=1&!is.na(Time_of_Day)&!(ID %in% collection.sets$ID)) %>% select(ID,Time_of_Day) %>% distinct() %>% group_by(ID) %>% count() %>% filter(n>1) %>% select(-n) %>% mutate(set="preanalytical: time of day") %>% bind_rows(collection.sets)
collection.sets <- d %>% filter(timespan.days<=1&!is.na(Sample_type)) %>% select(ID,Sample_type) %>% distinct()  %>% group_by(ID) %>% count() %>% filter(n>1) %>% select(-n) %>% mutate(set="preanalytical: tube type pairs") %>% bind_rows(collection.sets)
collection.sets <- d %>% filter(timespan.days<=1) %>% filter(ndraws>1&!ID %in% collection.sets$ID) %>% select(ID) %>% distinct() %>% mutate(set="preanalytical: other") %>% bind_rows(collection.sets)
collection.sets <- d %>% filter(timespan.days<=1&ndraws==1&!(ID %in% collection.sets$ID)) %>% select(ID) %>% distinct() %>% mutate(set="single sample") %>% bind_rows(collection.sets) 
collection.sets <- d %>% select(ID) %>% distinct() %>% full_join(collection.sets)
collection.sets <- d %>% filter(!is.na(replicate)) %>% select(ID) %>% distinct() %>% mutate(set="technical replicates") %>% bind_rows(collection.sets)

variable_input_dna <- d %>% filter(!is.na(ULP_input_DNA))  %>% select(ID,Plasma_ID,ULP_input_DNA) %>% distinct() %>% group_by(ID,Plasma_ID) %>% count() %>% filter(n>1) %>% inner_join(d %>% select(ID,Plasma_ID,ULP_input_DNA)) 
inputDNAstr <- paste("variable input DNA: ",min(variable_input_dna$ULP_input_DNA),"-",max(variable_input_dna$ULP_input_DNA)," ng",sep="")
collection.sets <- d %>% filter(ID %in% variable_input_dna$ID) %>% select(ID) %>% mutate(set=inputDNAstr) %>% bind_rows(collection.sets)
collection.sets <- collection.sets  %>% rename(collection.set=set) %>% ungroup()

tube.sets <- d %>% filter(Sample_type=="EDTA") %>% select(ID) %>% distinct() %>% mutate(tube.set="EDTA") 
tube.sets <- d %>% filter(Sample_type=="Streck") %>% select(ID) %>% distinct() %>% mutate(tube.set="Streck") %>% bind_rows(tube.sets) 

cbc.sets <- d %>% filter(CBC) %>% select(ID,Date_of_sample) %>% distinct() %>% group_by(ID) %>% count() %>% rename(ncbc=n) 
cbc.sets <- d  %>% select(ID) %>% distinct() %>% full_join(cbc.sets) %>% replace_na(list(ncbc=0))
cbc.sets <- cbc.sets %>% mutate(cbc.set=if_else(ncbc==0,"no CBC",if_else(ncbc==1,"CBC at one timepoint","CBC at more than one timepoint")))
cbc.sets <- cbc.sets %>% select(ID,cbc.set) %>% distinct() %>% ungroup()

# make data frame with meta data about dogs
dogs <- d %>% select(ID,Sex,Spay.neuter.status,Ancestry,Cancer_status,Cancer_type,Breed) %>% distinct()
# get date of first sample
dogs <- d %>% filter(!is.na(Date_of_sample)&!is.na(TF)&!is.na(DNA_conc)) %>% select(ID,Date_of_sample) %>% distinct() %>% group_by(ID) %>% summarise(Date_first_sample=min(Date_of_sample)) %>% right_join(dogs) 
# get age at first sample
dogs <- d %>% filter(!is.na(Age_at_sample)&!is.na(TF)&!is.na(DNA_conc)) %>% select(ID,Age_at_sample) %>% distinct() %>% group_by(ID) %>% summarise(Age_first_sample=min(Age_at_sample)) %>% right_join(dogs) 
dogs <- raw %>% mutate(P30=if_else(str_detect(Cohort,"P30"),TRUE,FALSE)) %>% select(ID,P30) %>% distinct() %>% right_join(dogs)


### Add first weight to dogs table (its complicated bc some dogs have no dates)
weight <- d %>% group_by(ID) %>% count() %>% filter(n==1) %>% inner_join(d) %>% select(ID,Date_of_sample,Weight_kg) %>% distinct()
weight <- d %>% filter(!(ID %in% weight$ID)&!is.na(Weight_kg)&!is.na(Date_of_sample)) %>% group_by(ID) %>% summarize(Date_of_sample=min(Date_of_sample)) %>% inner_join(d) %>% select(ID,Date_of_sample,Weight_kg)  %>% distinct() %>% bind_rows(weight) 
dogs <- dogs %>% left_join(weight %>% rename(Weight_first_sample=Weight_kg))
dogs <- dogs %>% select(-Date_of_sample) %>% distinct()

write.csv(dogs,"Table.Fig_2.csv",row.names=F)

## Make sequencing data type groups
raw <- raw %>% mutate(ULPWGS_data=if_else(ULPWGS_data_cf3|ULPWGS_data_cf4,TRUE,FALSE))
sequencing.sets <- raw %>% select(ID,Plasma_ID,ULPWGS_data,Normal_data,Tumor_data,WES_data) %>% distinct() %>% pivot_longer(c(ULPWGS_data,Normal_data,Tumor_data,WES_data)) %>% rename(fnd=value) %>% mutate(name=str_remove(name,"_data"))
sequencing.sets <- sequencing.sets %>% filter(fnd) %>% select(-Plasma_ID) %>% distinct()
sequencing.sets <- sequencing.sets %>% distinct() %>% arrange(ID,name) %>% pivot_wider(names_from=name,values_from=fnd) 
sequencing.sets <- sequencing.sets %>%  full_join(raw %>% select(ID) %>% distinct) %>% replace_na(list(ULPWGS=FALSE,WES=FALSE,Normal=FALSE,Tumor=FALSE))

sequencing.sets2 <- sequencing.sets %>% filter(ULPWGS) %>% mutate(set="Ultra lowpass (ULP)") %>% select(ID,set)
sequencing.sets2 <- sequencing.sets %>% filter(ULPWGS&Normal&Tumor) %>% mutate(set="ULP + paired tumor/normal(WGS)") %>% select(ID,set) %>% bind_rows(sequencing.sets2)
sequencing.sets2 <- sequencing.sets %>% filter(ULPWGS&Normal&WES&Tumor) %>% mutate(set="ULP + paired tumor/normal(WGS) + tumor (WES)") %>% select(ID,set) %>% bind_rows(sequencing.sets2)

sequencing.sets <- sequencing.sets2 %>% rename(sequencing.set=set) %>% mutate(facet="sequencing\ndata types")

highTF <- raw %>% filter(!is.na(TF)&TF>=0.1) %>% select(ID)  %>% distinct() %>% pull(ID)

sequencing.sets <- sequencing.sets %>% filter(sequencing.set=="ULP + paired tumor/normal(WGS)") %>% filter(ID %in% highTF) %>% mutate(sequencing.set=paste(sequencing.set,"(tumor fraction > 10%)")) %>% bind_rows(sequencing.sets)

cancer.sets <- d %>% select(ID,Cancer_type) %>% distinct() %>% group_by(Cancer_type) %>% count() %>% mutate(cancer.set=if_else(Cancer_type=="Healthy"|(n>3&Cancer_type!="Other"),Cancer_type,"Other")) %>% select(Cancer_type,cancer.set) %>% distinct()
cancer.sets <- d %>% select(ID,Cancer_type) %>% distinct() %>% full_join(cancer.sets) %>% select(-Cancer_type)

pd <- sequencing.sets %>% rename(set=sequencing.set) 
pd <- collection.sets %>% rename(set=collection.set) %>% mutate(facet="collection sets") %>% bind_rows(pd)
pd <- tube.sets %>% rename(set=tube.set) %>% mutate(facet="tube\ntypes") %>% bind_rows(pd)
pd <- cbc.sets %>% rename(set=cbc.set) %>% mutate(facet="CBC") %>% bind_rows(pd)

pd <- pd %>% left_join(cancer.sets) %>% select(ID,cancer.set) %>% rename(fill.set=cancer.set) %>% distinct() %>% full_join(pd)
pd <- pd %>% select(ID,fill.set,set,facet) %>% distinct() %>% group_by(fill.set,set,facet) %>% count() %>% rename(ndogs=n)

pdTotals <- pd %>% group_by(facet,set) %>% summarize(set.total.n=sum(ndogs))
levels <- pdTotals %>% group_by(facet,set) %>% summarize(set.total.n=max(set.total.n)) 

levels1 <- c("single sample","preanalytical: blood draw site","preanalytical: time of day","preanalytical: tube type pairs","preanalytical: other","longitudinal: short",longstr,"technical replicates",inputDNAstr,"EDTA","Streck","no CBC","CBC at one timepoint","CBC at more than one timepoint")
levels2 <- c("Ultra lowpass (ULP)","ULP + normal (WGS)","ULP + tumor (WGS)","ULP + tumor (WGS) (tumor fraction > 10%)","ULP + paired tumor/normal(WGS)","ULP + paired tumor/normal(WGS) (tumor fraction > 10%)","ULP + paired tumor/normal(WGS) + tumor (WES)")
levels <- tibble(set=c(rev(levels1),rev(levels2)))
levels <- levels %>% mutate(order=row_number())

fill.levels <- pd %>% filter(fill.set != "Healthy" & facet=="collection sets") %>% group_by(fill.set) %>% summarize(total=sum(ndogs))
fill.levels <- rev(c(fill.levels %>% arrange(total) %>% pull(fill.set),"Healthy"))

pd$set <- factor(pd$set,levels=levels$set)
pd$fill.set <- factor(pd$fill.set,levels=fill.levels)

pdTitles <- pd %>% filter(set==" ")
pd <- pd %>% filter(set != " ")
labels <- pd %>% ungroup() %>% select(facet,set) %>% distinct() %>% mutate(sumgrp=row_number())
labels <- pd %>% left_join(labels)
labels <- tibble(fill.set=rev(fill.levels),fill.order=c(1:length(fill.levels))) %>% right_join(labels)
labels <- labels %>% arrange(sumgrp,fill.order) %>% group_by(sumgrp) %>% mutate(xpos=cumsum(ndogs))
labels <- labels %>% ungroup() %>% rename(label=ndogs) %>% rename(ndogs=xpos) %>% select(set,fill.set,ndogs,label,facet)
labels <- labels %>% filter(label>=5) %>% mutate(xpos=((ndogs-label)+ndogs)/2)

samples_per_dog <- raw %>% select(ID,Plasma_ID) %>% distinct() %>% group_by(ID) %>% count() %>% full_join(collection.sets %>% select(ID,collection.set))
samples_per_dog <- samples_per_dog %>% ungroup() %>% group_by(collection.set) %>% summarize(min=min(n),max=max(n))
samples_per_dog <- samples_per_dog %>% mutate(range=if_else(min==max,if_else(min==1,paste(min,"sample/dog"),paste(min,"samples/dog")),paste(min,"-",max," samples/dog",sep=""))) 
samples_per_dog <- samples_per_dog %>% rename(set=collection.set) %>% select(set,range) 
##samples_per_dog <- samples_per_dog %>% mutate(range=if_else(set %in% c("preanalytical: EDTA tube","preanalytical: Streck tube"),"1 sample/dog",range))
ylabels <- pd %>% ungroup() %>% select(set) %>% distinct() %>% left_join(samples_per_dog)
ylabels <- ylabels %>% mutate(ylabel=if_else(is.na(range)|str_detect(set,"EDTA")|str_detect(set,"Streck"),set,paste(set," (",range,")",sep="")))

unique(pd$set)
p <- ggplot(pd,aes(x=ndogs,y=set)) + geom_bar(aes(fill=fill.set),width=0.6,stat="identity",position="stack")
p <- p + facet_grid(fct_relevel(facet,"sequencing\ndata types","collection sets","tube\ntypes","CBC")~.,scales="free_y",space="free") #,strip.position = "top")
p <- p + geom_text(aes(label=set.total.n,x=set.total.n+1),hjust=0,size=2.5,data=pdTotals)
p <- p + geom_text(aes(label=label,x=xpos),hjust=0.5,size=2,color="#FFFFFF",data=labels)
p <- p + scale_fill_manual(values=c("#525252","#d95f02","#1b9e77","#7570b3","#e7298a","#e6ab02","#a6761d","#66a61e")) ###d6604d","#fee0d2","#40004b","#9970ab","#c2a5cf","#e7d4e8","#1b7837","#a6dba0")))
p <- p + scale_x_continuous("# dogs",limits=c(0,155),breaks=c(0:10)*20)
p <- p + scale_y_discrete("",expand=c(0,0.5),breaks=ylabels$set,labels=ylabels$ylabel)
p <- p + theme_cowplot(12)
p <- p + theme(plot.title = element_blank(),
               axis.title.y=element_blank(),
               axis.title.x = element_text(size=9),axis.text.x = element_text(size=8),
               axis.text.y = element_text(size=7),legend.position = c(0.98,0.98),
               legend.text = element_text(size=6),legend.title = element_blank(),
               legend.key.size = unit(0.25, 'cm'),
               legend.justification.inside = c(1.05, 1.25),
               
               legend.box.just = "right",
               legend.margin = margin(2, 3, 3, 3),
               legend.background = element_rect(color="#969696"),
               strip.text.y = element_text(size=6,hjust=0.5,face="bold",color="#f7f7f7"),
               strip.background=element_rect(fill="#525252"))
p_nsets_samplesets <- p

#### Make second panel - empty for now 

#p_concordance <- ggplot() + theme_void()

##### Make third panel - overlap


raw <- as_tibble(read.csv(infile_conmuts,header=T,sep="\t",na.strings=c("","NA","NI","UNK"),strip.white=TRUE,stringsAsFactors=TRUE)) %>% filter(FILTER=="PASS") %>% select(-FILTER)
raw <- raw %>% select(ID,Sample_type,Sample_ID,CHROM,POS,Var_ID,Clonal)
locations <- raw %>% select(CHROM,POS) %>% distinct() %>% arrange(CHROM,POS) %>% mutate(mutID=row_number())
raw <- raw %>% inner_join(locations)

# Plot overlap of cfDNA and tumor DNA
d1 <- raw %>% select(-Sample_ID,-Clonal) %>% mutate(value=TRUE) %>% pivot_wider(names_from=Sample_type) %>% replace_na(list(tumor=FALSE,cfDNA=FALSE))
d2 <- raw %>% select(-Sample_ID) %>% pivot_wider(names_from=Sample_type,values_from=Clonal)  %>% replace_na(list(tumor=FALSE,cfDNA=FALSE)) %>% rename(tumor.cl=tumor,cfDNA.cl=cfDNA)
d <- d1 %>% full_join(d2)

order_match <- tibble(matchset=c("both","cfDNA","tumor","neither"),morder=c(1:4)*10)
order_clonal <- tibble(clonalset=c("clonal","subclonal","subclonal in cfDNA","subclonal in tumor"),corder=c(1:4))

d <- d %>% mutate(matchset=if_else(tumor&cfDNA,"both",if_else(tumor&!cfDNA,"tumor",if_else(cfDNA&!tumor,"cfDNA","neither"))))
d <- d %>% mutate(clonalset=if_else(tumor&!cfDNA|cfDNA&!tumor,if_else(tumor.cl|cfDNA.cl,"clonal","subclonal"),if_else(tumor&cfDNA,if_else(tumor.cl&cfDNA.cl,"clonal",if_else(!tumor.cl&!cfDNA.cl,"subclonal",if_else(tumor.cl&!cfDNA.cl,"subclonal in tumor",if_else(!tumor.cl&cfDNA.cl,"subclonal in cfDNA","error")))),"error2")))
d <- d %>% left_join(order_match) %>% left_join(order_clonal) %>% mutate(order=morder+corder)
d <- d %>% select(-tumor,-cfDNA,-tumor.cl,-cfDNA.cl,-morder,-corder)
d <- d %>% mutate(baseset=matchset,matchset=paste(matchset," (",clonalset,")",sep="")) %>% select(-clonalset)

mlevels <- d %>% select(matchset,order) %>% distinct() %>% arrange(-order) %>% pull(matchset)

pd <- d %>% group_by(ID,matchset,baseset) %>% count()
totals <- pd %>% group_by(ID) %>% summarize(total=sum(n)) 
pd <- totals %>% inner_join(pd) %>% mutate(frac=n/total)

pd <- pd %>% mutate(label=paste(round(frac*100,0),"%",sep=""))
pd <- pd %>% arrange(ID,desc(matchset)) %>% group_by(ID) %>% mutate(sum=cumsum(n),row=row_number())

totals <- totals %>% mutate(string=paste("N=",total,sep=""))
starts <- pd %>% select(ID,sum,row) %>% mutate(row=row+1) %>% rename(start=sum)
pd <- pd %>% left_join(starts) %>% replace_na(list(start=0))
pd <- pd %>% mutate(labelx=(start+sum)/2) 
pd <- pd %>% select(ID,baseset,matchset,n,total,frac,label,labelx)

ylevels <- pd %>% select(ID,total) %>% distinct() %>% arrange(total) %>% pull(ID)

pd$matchset <- factor(pd$matchset,levels=mlevels)
pd$ID <- factor(pd$ID,levels=ylevels)

subtotals <- pd %>% group_by(ID,baseset) %>% summarize(subtotal=sum(n))
subtotals <- subtotals %>% group_by(ID) %>% mutate(xend=cumsum(subtotal))
subtotals <- subtotals %>% mutate(xpos=(xend-subtotal+xend)/2)
subtotals <- subtotals %>% left_join(totals %>% select(-string)) %>% mutate(percent=paste(round(subtotal/total*100,0),"%",sep="")) %>% mutate(matchset=paste(baseset,"clonal"))

p <- ggplot(pd,aes(x=n,y=ID)) + geom_bar(aes(fill=matchset),stat="identity",width=0.5)
p <- p + geom_bar(aes(x=total),color="#1a1a1a",fill=NA,data=totals,stat="identity",width=0.5,linewidth=0.5)
p <- p + geom_text(aes(x=xpos,label=percent),color="#67001f",hjust=0.5,vjust=0,nudge_y = 0.35,size=2,data=subtotals %>% filter(baseset=="both"))
p <- p + geom_text(aes(x=xpos,label=percent),color="#2166ac",hjust=0.5,vjust=0,nudge_y = 0.35,size=2,data=subtotals %>% filter(baseset=="cfDNA"))
p <- p + geom_text(aes(x=xpos,label=percent),color="#4d4d4d",hjust=0.5,vjust=0,nudge_y = 0.35,size=2,data=subtotals %>% filter(baseset=="tumor"))
p <- p + geom_text(aes(x=total,label=string),hjust=0,vjust=0.5,data=totals,nudge_x=0.6,size=2)
p <- p + scale_y_discrete("") #,breaks=ylabels$ID,labels=ylabels$ID)
p <- p + scale_x_continuous("# somatic mutations",,expand = expansion(mult = c(0.05,0.2))) #,limits=c(0,120),
p <- p + scale_fill_manual(values=c("#bababa","#4d4d4d","#4393c3","#2166ac","#f4a582","#d6604d","#b2182b","#67001f"))
p <- p + scale_color_manual(values=c("#bababa","#4d4d4d","#4393c3","#2166ac","#f4a582","#d6604d","#b2182b","#67001f"))
p <- p + theme_cowplot(12)
p <- p + theme(plot.title = element_text(hjust=0,size=7,face="bold"),
               plot.subtitle = element_text(hjust=0,size=7),
               axis.line.y=element_blank(),axis.ticks.y=element_blank(),
               legend.title = element_blank(),
               axis.title = element_text(size=7),axis.text = element_text(size=6),
              legend.key.spacing.y=unit(0.1, 'cm'),
               legend.position="bottom", #c(0.98,0.02),
               legend.justification = c("right", "bottom"),
               legend.text = element_text(size=4.5,vjust=0.5),legend.key.size = unit(0.2, 'cm'),
               strip.background=element_rect(fill="#FFFFFF"),
               strip.text = element_text(size=7,hjust=0.5,face="bold",color="#252525"),
               legend.box.background = element_rect(color="#d9d9d9", size=0.5),legend.margin = margin(3,5,5,5))

p_overlap <- p

pd <- as_tibble(read.csv(infile_depth,header=T,sep="\t",na.strings=c("","NA","NI","UNK"),strip.white=TRUE,stringsAsFactors=TRUE)) 
pd <- pd %>% filter(!is.na(tumor_log2_ratio)&!is.na(plasma_log2_ratio))
pd <- pd %>% inner_join(samples_over_10per)
max_xy <- max(pd$tumor_log2_ratio,pd$plasma_log2_ratio)
min_xy <- min(pd$tumor_log2_ratio,pd$plasma_log2_ratio)
labels <- c(0.125,0.25,0.5,1,2,4,8,16)

p <- ggplot(pd,aes(x=tumor_log2_ratio,y=plasma_log2_ratio))
p <- p + rasterise(geom_point(color="#2166ac",size=0.5,alpha=0.5,shape=16),dpi=600)
p <- p + geom_abline(color="#4d4d4d",alpha=0.75)
p <- p + scale_x_continuous("read depth (tumor)",limits=c(min_xy,max_xy),labels=labels,breaks=log2(labels))
p <- p + scale_y_continuous("read depth (plasma)",limits=c(min_xy,max_xy),labels=labels,breaks=log2(labels))
p <- p + stat_cor(
  aes(label = sprintf("R = %0.2f, p = %g", after_stat(r), after_stat(p))), 
  method = "spearman", 
  size = 3,
  output.type = "text"
)
p <- p + theme_cowplot(12)
p <- p + theme(axis.title = element_text(size=7),axis.text = element_text(size=6),
               axis.line=element_line(linewidth=0.25),axis.ticks=element_line(linewidth=0.25),
               legend.position="none")
p_depth <- p

row1 <- plot_grid(p_nsets_samplesets,ncol=1,label_size = 12,labels = LETTERS[1])
row2 <- plot_grid(p_depth,p_overlap,ncol=2,label_size = 12,labels = LETTERS[-1],rel_widths = c(0.5,0.6))
grid <- plot_grid(row1,row2,ncol=1,rel_heights = c(1.5,1))
ggsave(grid,filename=pdfname,width=6.5,height=6)

