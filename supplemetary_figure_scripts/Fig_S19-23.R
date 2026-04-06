library(tidyverse)
library(cowplot)
library(ggpubr)
library(rstatix)
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

pdfname <- paste(pdfdir,"Fig_SLONGITUDINAL.pdf",sep="")

rename_stats <- tibble(order=c(2,1,3),statistic=c("Tumor fraction","cfDNA (ng/mL)","Fragment size ratio"),stat=c("TF","log_DNA_conc","Fragment_size_ratio") ) %>% arrange(order) %>% select(-order)
cohorts <- c("OS_Longitudinal","LSA_U01_treatment_cohort_2","LSA_U01_treatment_cohort_1")

raw <- as_tibble(read.csv(infile,header=T,sep="\t",na.strings = c("NA","N/A","","UNK")))
raw <- raw %>% filter(Cohort %in% cohorts)

raw <- raw %>% mutate(Date_of_sample=str_replace(Date_of_sample," \\(approx\\)","")) %>% mutate(Date_of_sample=parse_date_time(Date_of_sample,"m/d/y"))
raw <- raw %>% mutate(Date_of_diagnosis=str_remove(Date_of_diagnosis,"10/28/2011 and ")) %>% mutate(Date_of_diagnosis=str_replace(Date_of_diagnosis," \\(approx\\)","")) %>% mutate(Date_of_diagnosis=parse_date_time(Date_of_diagnosis,"m/d/y"))
raw <- raw %>% mutate(Date_of_progressive_disease2=parse_date_time(Date_of_progressive_disease,"m/d/y"),Last_known_date_if_still_alive=parse_date_time(Last_known_date_if_still_alive,"m/d/y"))
raw <- raw %>% mutate(log_DNA_conc=log10(DNA_conc)) #%>% select(-DNA_conc)


# convert Disease status in readable form and align terminology for osteosarcoma and lymphoma
Disease_status <- tibble(  Disease_status=c("disease present","PR","NED","PD"),  Disease_status_new=c("disease present","partial resp.","no disease","progressive disease"))
Disease_status <- raw %>% select(Disease_status) %>% distinct() %>% inner_join(Disease_status)
Disease_status$status_order <- c(4,3,2,1)
raw <- raw %>% left_join(Disease_status) %>% rename(Disease_status_short=Disease_status) %>% rename(Disease_status=Disease_status_new)
Disease_status <- raw %>% select(Disease_status,status_order) %>% distinct()

## make everything references on Date of first sample
raw <- raw %>% group_by(ID) %>% summarize(Date_of_first_sample=min(Date_of_sample)) %>% full_join(raw)
raw <- raw %>% mutate(timepoint=as.Date(Date_of_sample)-as.Date(Date_of_first_sample))
raw <- raw %>% mutate(OST=OST-(as.numeric(as.Date(Date_of_first_sample)-as.Date(Date_of_diagnosis))))
raw <- raw %>% mutate(PFS=PFS-(as.numeric(as.Date(Date_of_first_sample)-as.Date(Date_of_diagnosis))))

### keep the samples with the largest input DNA amount
inputdna <- raw %>% filter(!is.na(ULP_input_DNA)) %>% select(ID,Sample_ID,Plasma_ID,ULP_input_DNA) %>% group_by(ID,Plasma_ID) %>% summarize(max=max(ULP_input_DNA)) %>% inner_join(raw,relationship = "many-to-many")
inputdna <- inputdna %>% ungroup() %>%  mutate(keep=if_else(ULP_input_DNA==max,TRUE,FALSE)) %>% select(ID,Sample_ID,keep)
raw <- raw %>% left_join(inputdna,relationship = "many-to-many") %>% filter(is.na(ULP_input_DNA)|keep) %>% select(-keep)

### separate out dog level metrics
dogs <- raw %>% select(ID,Cancer_type,ID,Date_of_birth, Date_of_diagnosis, Date_of_progressive_disease, PFS,Date_of_death, Last_known_date_if_still_alive, OST, Sex, Spay.neuter.status, Breed,Primary_diagnosis,Cohort,Institution) %>% distinct()
dogs <- dogs %>% mutate(longname=if_else(is.na(ID)|ID==ID,ID,paste(ID,"\n",ID,sep="")))

raw <- raw %>% select(ID,Sample_ID,Plasma_ID,replicate,ULP_input_DNA,log_DNA_conc,TF,Fragment_size_ratio,timepoint,Date_of_sample,Age_at_sample,Weight_kg,Cancer_type,Disease_status) %>% distinct()
### add sample order number colun
raw <- raw %>% mutate(days=as.numeric(timepoint))
### mark first instance of progressive disease
progressive <- raw %>% filter(Disease_status=="progressive disease") %>% select(ID,timepoint) %>% group_by(ID) %>% summarize(timepoint=min(timepoint)) %>% mutate(relapse=TRUE)
raw <- raw %>% left_join(progressive) %>% replace_na(list(relapse=FALSE))

ost_dfi <- dogs %>% select(ID,PFS,OST) %>% pivot_longer(c(PFS,OST)) %>% filter(!is.na(value)) %>% rename(stat=name,days=value) %>% mutate(value=0)

raw_order <- raw %>% select(ID,days) %>% bind_rows(ost_dfi %>% select(ID,days)) %>% distinct()
raw_order <- raw_order %>% arrange(ID,days) %>% group_by(ID) %>% mutate(order=row_number())

raw <- raw %>% left_join(raw_order)
ost_dfi <- ost_dfi %>% left_join(raw_order)

#### Make boxplot comparing metrics by disease status
colors4 <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3")


pdAll <- raw %>% select(ID,Sample_ID,Plasma_ID,timepoint,order,log_DNA_conc,TF,Fragment_size_ratio,Disease_status,relapse) %>% pivot_longer(c(log_DNA_conc,TF,Fragment_size_ratio))  %>% rename(stat=name)
pdAll <- dogs %>% select(ID,Cancer_type,longname,PFS,OST) %>% inner_join(pdAll)

pdAll <- pdAll %>% mutate(plotbin=if_else(timepoint==0,"disease present",Disease_status))
                                          #if_else(Disease_status %in% c("disease present","progressive disease","partial resp."),Disease_status,
                                           #       paste(Disease_status,if_else(timepoint<30,"first_month","later"),sep=".")))) %>% mutate(plotbin=str_replace(plotbin,"\\.\\.","."))
pdAll <- pdAll %>% left_join(rename_stats)

## remove partial response because we don't have it for osteosaroma 
pdAll <- pdAll %>% filter(Disease_status!="partial resp.")

pd_progressive <- pdAll %>% filter(relapse)

#plotbins <- nextUp %>% select(Sample_ID,plotbin)
plotbins <- c("disease present","no disease","progressive disease")
comparisons <- list(c(plotbins[1],plotbins[2]),c(plotbins[2],plotbins[3])) #,c(plotbins[3],plotbins[4]))

pdAll$plotbin <- factor(pdAll$plotbin,levels=plotbins)
xbreaks <- plotbins
xlabels <- c("disease\npresent","no\ndisease","progressive\ndisease")

theme_boxes <- function(){ 
  theme_cowplot(12) %+replace%
    theme(     plot.title = element_text(hjust=0,size=7,face="bold",margin=margin(0,0,4,0)),
               plot.subtitle = element_text(hjust=0,size=6,margin=margin(0,0,4,0),lineheight = 1),
               axis.title.x=element_blank(),
               axis.title.y = element_text(size=8,angle=90,vjust=1),
               axis.text.x = element_text(size=4.5,vjust=1,margin=margin(2,0,0,4),hjust=0.5),
               axis.text.y = element_text(size=5,margin=margin(2,0,0,4)),
               strip.text = element_text(size=6,hjust=0.5,face="bold",color="#f7f7f7"),
               strip.background=element_rect(fill="#525252"),
               panel.border = element_rect(color = "#525252", fill = NA, size = 0.5),
               legend.position="none")
}

plot_box <- function(statIn) {
  
pd <- pdAll %>% filter(stat==statIn) %>% filter(!is.na(value)) 

p <- ggplot(pd,aes(x=plotbin,y=value))
p <- p + geom_boxplot(color="#525252",outlier.size=1,outlier.shape=NA)
p <- p + geom_jitter(aes(color=Cancer_type),shape=16,size=1,alpha=0.5,height=0,width=0.2)
p <- p + stat_compare_means(method = "t.test",comparisons= comparisons,size=1.5)
p <- p + scale_color_manual(values=c("#ef3b2c","#8073ac"))
p <- p + facet_wrap(~Cancer_type,nrow=1,scales="free")
p <- p + scale_x_discrete("",breaks=xbreaks,labels=xlabels)
if (statIn=="log_DNA_conc") {
  p <- p + scale_y_continuous(unique(pd$statistic),breaks=c(0:4),labels=10**c(0:4)) 
} else {
  p <- p + scale_y_continuous(unique(pd$statistic))
}
p <- p + ggtitle(unique(pd$statistic),subtitle="t-test")
p <- p + theme_boxes()
p
}
p_conc <- plot_box("log_DNA_conc")
p_TF <- plot_box("TF")
p_frag <- plot_box("Fragment_size_ratio")

pdMean <- pdAll %>% group_by(ID,Cancer_type,plotbin,stat) %>% summarize(value=mean(value)) %>% left_join(rename_stats)
pdMean$plotbin <- factor(pdMean$plotbin,levels=plotbins) 

plot_paired <- function(statIn) {
  statIn <- "TF"
pd <- pdMean %>% filter(stat==statIn) %>% filter(!is.na(value)) 
compare_set1 <- pd %>% filter(plotbin %in% plotbins[1:2]) %>% group_by(ID) %>% count() %>% filter(n==2) %>% select(ID) 
compare_set2 <- pd %>% filter(plotbin %in% plotbins[2:3]) %>% group_by(ID) %>% count() %>% filter(n==2) %>% select(ID) 
okIds <- compare_set1 %>% bind_rows(compare_set2) %>% distinct()
pd <- pd %>% inner_join(okIds)
pd <- pd %>% group_by(ID) %>% count() %>% filter(n>=2) %>% select(ID) %>% inner_join(pd)
}

p <- ggplot(pd,aes(x=plotbin,y=value))
p <- p + geom_point(aes(color=Cancer_type),shape=16,size=1,alpha=0.9)
p <- p + geom_line(aes(group=ID),color="#707070",alpha=0.5,size=0.5)
p <- p + stat_compare_means(method = "t.test",paired=TRUE,size=1.5,comparisons=list(plotbins[1:2]),data=pd %>% inner_join(compare_set1)) #,data=pd %>% filter(ID %in% compare_set1))
p <- p + stat_compare_means(method = "t.test",paired=TRUE,size=1.5,comparisons=list(plotbins[2:3]),data=pd %>% inner_join(compare_set2)) #,data=pd %>% filter(ID %in% compare_set1))
p <- p + scale_color_manual(values=c("#ef3b2c","#8073ac"))
p <- p + facet_wrap(~Cancer_type,nrow=1,scales="free")
p <- p + scale_x_discrete("",breaks=xbreaks,labels=xlabels)
if (statIn=="log_DNA_conc") {
  p <- p + scale_y_continuous(unique(pd$statistic),breaks=c(0:4),labels=10**c(0:4)) 
} else {
  p <- p + scale_y_continuous(unique(pd$statistic))
}
p <- p + ggtitle(unique(pd$statistic),subtitle="paired t-test")
p <- p + theme_boxes()
p


p_conc2 <- plot_paired("log_DNA_conc")
p_TF2 <- plot_paired("TF")
p_frag2 <- plot_paired("Fragment_size_ratio")


row1 <- plot_grid(p_conc,p_TF,p_frag,nrow=1,labels=LETTERS[1:3],rel_widths = c(1,2,2))
row2 <- plot_grid(p_conc2,p_TF2,p_frag2,nrow=1,labels=LETTERS[4:6],rel_widths = c(1,2,2))
grid <- plot_grid(row1,row2,ncol=1)
ggsave(grid,filename=pdfname,width=6.5,height=5)


theme_long <- function(){ 
  theme_cowplot(12) %+replace%
    theme(     plot.title = element_text(hjust=0,size=8,face="bold",margin=margin(0,0,4,0)),
               plot.subtitle = element_text(hjust=0,size=7,margin=margin(0,0,4,0),lineheight = 1),
               axis.title.x=element_blank(),
               axis.line.x=element_blank(),
               axis.line.y=element_line(color="#878787"),
               axis.title.y = element_text(size=8,angle=90,vjust=1),
               axis.text.x = element_blank(),
               axis.text.y = element_text(size=5), 
               axis.ticks.x = element_blank(),
               strip.text = element_text(size=6,hjust=0.5,face="bold",color="#f7f7f7",lineheight=1.1),
               strip.background=element_rect(fill="#525252"),
               panel.border = element_blank(),
               legend.position="bottom",
               #legend.justification = c("right", "bottom"),
               legend.margin = margin(1,1,1,1),
               legend.text = element_text(size=4.5,vjust=0.5),legend.key.size = unit(0.20, 'cm'),
               legend.title = element_blank(),
               legend.box.background = element_rect(color="white", size=0.5))
  
}

statuses <- tibble(Disease_status=c("no disease","partial resp.","disease present","progressive disease"),
                       status_color=c("#4d4d4d","#fc9272","#ef3b2c","#a50f15"),
                       status_order=c(1:4))


pdLongAll <- raw %>% select(ID,Sample_ID,Plasma_ID,timepoint,order,log_DNA_conc,TF,Fragment_size_ratio,Disease_status,relapse) %>% pivot_longer(c(log_DNA_conc,TF,Fragment_size_ratio))  %>% rename(stat=name) %>% mutate(days=as.numeric(timepoint))
pdLongAll <- dogs %>% select(ID,Cancer_type,longname) %>% inner_join(pdLongAll) %>% left_join(rename_stats)

pdPFSOST_all <- dogs %>% select(ID,Cancer_type,longname) %>% inner_join(ost_dfi) 
relpos <- pdLongAll %>% select(ID,days,order) %>% distinct()
pdLongAll <- pdLongAll %>% left_join(relpos %>% rename(PFS=days,PFS_order=order),relationship = "many-to-many")
pdLongAll <- pdLongAll %>% left_join(relpos %>% rename(OST=days,OST_order=order),relationship = "many-to-many")

#statIn <- "TF"
#cancerTypeIn <- "Osteosarcoma"

plot_long <- function(statIn,cancerTypeIn) {
  
pd <- pdLongAll %>% filter(stat==statIn&Cancer_type==cancerTypeIn)

pd_status <- pd %>% ungroup() %>% select(Disease_status) %>% distinct() %>% inner_join(statuses) %>% arrange(status_order)
pd$Disease_status <- factor(pd$Disease_status,levels=pd_status$Disease_status)


pdPFSOST <- pdPFSOST_all %>% filter(Cancer_type==cancerTypeIn)
pd_w_allpts <- pd %>% bind_rows(pdPFSOST)
labels <- pd_w_allpts %>% ungroup() %>% select(order,value,longname,days,order) %>% distinct()

ncol <- 5 
ntot <- length(unique(pd$ID))
nrow <- ceiling(ntot/ncol)
miny = min(pd$value)
maxy = max(pd$value)

p <- ggplot(pd_w_allpts,aes(x=order,y=value))
p <- p + geom_hline(yintercept = miny,linewidth=0.5,color="#878787")
p <- p + geom_segment(aes(x=order,xend=order),y=miny,yend=maxy,size=2.5,alpha=0.25,color="#a50f15",data=pdPFSOST %>% filter(stat=="PFS"))
p <- p + geom_segment(aes(x=order,xend=order),y=miny,yend=maxy,size=0.5,color="#4d4d4d",data=pdPFSOST %>% filter(stat=="OST"),linetype = 2)
p <- p + geom_line(color="#969696",linewidth=1,alpha=0.5,data=pd)
p <- p + geom_point(aes(color=Disease_status),size=1.5,data=pd)
p <- p + geom_point(aes(color=Disease_status),shape=21,size=3.5,data=pd %>% filter(relapse)) #,shape=plot_shape_bin),size=3.5,data=pd %>% filter(plot_shape_bin!="other"))
p <- p + geom_text(aes(label=days),y=miny-((maxy-miny)/10),hjust=0.5,vjust=1,size=1.75,data=labels)
p <- p + facet_wrap(~longname,ncol=ncol,scales="free_x")
p <- p + ggtitle(unique(pd$statistic))
if (statIn=="log_DNA_conc") {
  p <- p + scale_y_continuous(unique(pd$statistic),breaks=c(0:4),labels=10**c(0:4),expand = expansion(mult = c(0.25, .25)))
} else {
  p <- p + scale_y_continuous(unique(pd$statistic),expand = expansion(mult = c(0.25, .25)))
}

p <- p + scale_x_continuous("",expand = expansion(mult = c(0.1,0.1)))
p <- p + scale_color_manual(values=pd_status$status_color)
#p <- p + scale_fill_manual("dog",values=colors4)
p <- p + scale_shape_manual("dog",values=c(21,23))
p <- p + theme_long()
pdfout <- paste(pdfdir,"Fig_S19.",cancerTypeIn,".",statIn,".pdf",sep="")
ggsave(p,filename=pdfout,width=6.5,height=nrow)

}

plot_long("log_DNA_conc","Lymphoma")
plot_long("TF","Lymphoma")
plot_long("TF","Osteosarcoma")
plot_long("Fragment_size_ratio","Lymphoma")
plot_long("Fragment_size_ratio","Osteosarcoma") 


