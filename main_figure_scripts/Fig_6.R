library(tidyverse)
library(cowplot)
# This package automatically finds the root folder of the downloaded github project.
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ==============================================================================
# 1. SETUP & FILE PATHS
# ==============================================================================

# Define the input file
infile <- here("data", "BB_Supplementary_Data_1.txt")

# Define the output directory
pdfdir <- here("figures")

# Create the output directory if it doesn't already exist
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir, recursive = TRUE)
  message(paste("Created output directory:", pdfdir))
}

rename_stats <- tibble(order=c(2,1,3),statName=c("Tumor fraction","cfDNA concentration","Fragment size ratio"),name=c("TF","log_DNA_conc","Fragment_size_ratio") ) %>% arrange(order)

raw <- as_tibble(read.csv(infile,header=T,sep="\t",na.strings = c("NA","N/A","","UNK","?")))
raw$X <- NULL
raw <- raw %>% filter(Cohort=="ARM")

d <- raw %>% select(ID,Sample_ID,Plasma_ID,First_Sample,DNA_conc,TF,Fragment_size_ratio,Sample_time_point,Tumor_size_cm,replicate) %>% ungroup()

timepoints <- d %>% select(Sample_time_point) %>% distinct() %>% mutate(day=str_remove(str_remove(Sample_time_point," Hour [0-9]+"),"Day "))
timepoints <- timepoints %>% mutate(hour=str_remove(str_remove(Sample_time_point,"Day [0-9]+ "),"Hour "))
timepoints <- timepoints %>% mutate(hour=as.numeric(hour),day=as.numeric(day))
timepoints <- timepoints %>% arrange(day+(hour/24)) %>% mutate(timepoint=row_number()-1,label=str_replace(Sample_time_point,"  ","\n"))
timepoints <- timepoints %>% mutate(label=str_replace(label," Hour","\nHr"))
timepoints <- timepoints %>% mutate(label=str_replace(label,"Day","D"))

d <- d %>% left_join(timepoints)

d <- d %>% mutate(log_DNA_conc=log10(DNA_conc)) %>% select(-DNA_conc)
d <- d %>% pivot_longer(c(log_DNA_conc,TF,Fragment_size_ratio)) %>% filter(!is.na(value)) %>% rename(stat=name)

pdAll <- d %>% select(ID,timepoint,stat,value)  %>% left_join(rename_stats) %>% distinct()
replicate_number <- pdAll %>% ungroup() %>% select(ID,timepoint,stat,value) %>% group_by(ID,timepoint,stat) %>% mutate(replicate=row_number())
pdAll <- pdAll %>% full_join(replicate_number) %>% filter(replicate==1)

#### Make lineplots comparing metrics by timepoint
theme_short <- function(){ 
  theme_cowplot(12) %+replace%
    theme(     plot.title = element_text(hjust=0,size=8,face="bold",margin=margin(0,0,4,0)),
               plot.subtitle = element_text(hjust=0,size=7,margin=margin(0,0,4,0),lineheight = 1),
               axis.title.x=element_blank(),
               axis.title.y = element_text(size=7,angle=90,vjust=1,face="bold",margin=margin(0,4,0,0)),
               axis.text.x = element_text(size=6),
               axis.text.y = element_text(size=6),
               legend.position="none")
}
#colors3 <- c("#66c2a5","#8da0cb","#e78ac3","#fc8d62")
colors3 <- c("#66c2a5","#e78ac3","#8da0cb") #"#fc8d62")
fills4 <- c("#66c2a5","#e78ac3","#8da0cb","#fc8d62")
shapes4 <- c(22,23,24,21)
ids_ordered <- c("ARM_412703","ARM_409177","ARM_S226485","ARM_408932")
response_ordered <- c("poor","good","partial")

response_type <- tibble(ID=c("ARM_412703","ARM_S226485","ARM_408932","ARM_409177"),response=c("poor","partial","partial","good"))
pdAll <- pdAll %>% left_join(response_type)
### Make Concentration
pd <- pdAll %>% filter(stat=="log_DNA_conc") 
pd$ID <- factor(pd$ID,levels=ids_ordered)
pd$response <- factor(pd$response,levels=response_ordered)

doglabels <- pd %>% group_by(stat,ID) %>% summarize(timepoint=max(timepoint)) %>% inner_join(pd)
title <- paste("cfDNA concentration\n",sep=" ")

##### Change relative timepoint 
relative_timepoint <- 0
timepoint_description <- "first timepoint"

pdRel <- pdAll %>% ungroup() %>% filter(timepoint==relative_timepoint) %>% select(ID,stat,value) %>% rename(firstvalue=value)

pdRelAll <- pdAll %>% left_join(pdRel) %>% mutate(value=value/firstvalue) %>% select(-firstvalue)

### Make Relative Concentration
pd <- pdRelAll %>% filter(stat=="log_DNA_conc") 

pd$ID <- factor(pd$ID,levels=ids_ordered)
pd$response <- factor(pd$response,levels=response_ordered)

doglabels <- pd %>% group_by(stat,ID) %>% summarize(timepoint=max(timepoint)) %>% inner_join(pd)
title <- paste("cfDNA concentration change\n(relative to ",timepoint_description,")",sep="")

p <- ggplot(pd,aes(x=timepoint,y=value))
p <- p + geom_vline(xintercept =c(3),color="#fb8072",linetype=1,size=3,alpha=0.3)

p <- p + geom_vline(xintercept =c(2.5,6.5),color="grey30",linetype=2)
p <- p + geom_hline(yintercept=1,color="#535353",linewidth=0.25)

p <- p + geom_line(aes(group=ID,color=response),linewidth=0.5)
p <- p + geom_text(aes(label=response,color=response),hjust=0,data=doglabels,nudge_x=0.2,size=2)
p <- p + ggtitle(title)
p <- p + geom_point(aes(fill=ID,shape=ID,y=value,color=response),alpha=0.75)
p <- p + scale_y_continuous(str_replace(unique(pd$statistic),"concentration","conc."))
p <- p + scale_x_continuous("",breaks=timepoints$timepoint,label=timepoints$label,,expand = expansion(mult = c(0.05,0.1)))
p <- p + scale_color_manual("dog",values=colors3)
p <- p + scale_fill_manual("dog",values=fills4)
p <- p + scale_shape_manual("dog",values=shapes4)

p <- p + theme_short()
p_rel_conc <- p


### Make Relative TF
pd <- pdRelAll %>% filter(stat=="TF") 
pd$ID <- factor(pd$ID,levels=ids_ordered)
pd$response <- factor(pd$response,levels=response_ordered)

doglabels <- pd %>% group_by(stat,ID) %>% summarize(timepoint=max(timepoint)) %>% inner_join(pd)
title <- paste("Tumor fraction change\n(relative to ",timepoint_description,")",sep="")

p <- ggplot(pd,aes(x=timepoint,y=value))
p <- p + geom_vline(xintercept=c(3),color="#fb8072",linetype=1,size=3,alpha=0.3)

p <- p + geom_vline(xintercept =c(2.5,6.5),color="grey30",linetype=2)
p <- p + geom_hline(yintercept=1,color="#535353",linewidth=0.25)

p <- p + geom_line(aes(group=ID,color=response),linewidth=0.5)
p <- p + geom_text(aes(label=response,color=response),hjust=0,data=doglabels,nudge_x=0.2,size=2)
p <- p + ggtitle(title)
p <- p + geom_point(aes(fill=ID,shape=ID,y=value,color=response),alpha=0.75)
p <- p + scale_y_continuous(unique(pd$statistic))
p <- p + scale_x_continuous("",breaks=timepoints$timepoint,label=timepoints$label,expand = expansion(mult = c(0.05,0.1)))
p <- p + scale_color_manual("dog",values=colors3)
p <- p + scale_fill_manual("dog",values=fills4)
p <- p + scale_shape_manual("dog",values=shapes4)
p <- p + theme_short()
p_rel_TF <- p

### Make Relative Fragment size ratio
pd <- pdRelAll %>% filter(stat=="Fragment_size_ratio") 
pd$ID <- factor(pd$ID,levels=ids_ordered)
pd$response <- factor(pd$response,levels=response_ordered)

doglabels <- pd %>% group_by(stat,ID) %>% summarize(timepoint=max(timepoint)) %>% inner_join(pd)
title <- paste("Fragment size ratio change\n(relative to ",timepoint_description,")",sep="")

p <- ggplot(pd,aes(x=timepoint,y=value))
p <- p + geom_vline(xintercept =c(3),color="#fb8072",linetype=1,size=3,alpha=0.3)
p <- p + geom_vline(xintercept =c(2.5,6.5),color="grey30",linetype=2)

p <- p + geom_hline(yintercept=1,color="#535353",linewidth=0.25)
p <- p + geom_line(aes(group=ID,color=response),linewidth=0.5)
p <- p + geom_text(aes(label=response,color=response),hjust=0,data=doglabels,nudge_x=0.2,size=2)
p <- p + ggtitle(title)
p <- p + geom_point(aes(fill=ID,shape=ID,y=value,color=response),alpha=0.75)
p <- p + scale_y_continuous(unique(pd$statistic))
p <- p + scale_x_continuous("",breaks=timepoints$timepoint,label=timepoints$label,,expand = expansion(mult = c(0.05,0.1)))
p <- p + scale_color_manual("dog",values=colors3)
p <- p + scale_fill_manual("dog",values=fills4)
p <- p + scale_shape_manual("dog",values=shapes4)
p <- p + theme_short()
p_rel_frag <- p

### Make Lymphnodes
pd <- d  %>% select(ID,timepoint,Tumor_size_cm) %>% distinct() %>% filter(!is.na(Tumor_size_cm))
pd <- pd %>% left_join(response_type)
pd <- pd %>% filter(timepoint==0) %>% rename(rel_size=Tumor_size_cm) %>% select(-timepoint) %>% full_join(pd) %>% mutate(rel_size=Tumor_size_cm-rel_size)
doglabels <- pd %>% filter(!is.na(Tumor_size_cm)) %>% group_by(ID) %>% summarize(timepoint=max(timepoint)) %>% inner_join(pd)
pd$ID <- factor(pd$ID,levels=ids_ordered)
pd$response <- factor(pd$response,levels=response_ordered)

title <- paste("Lymph node size (relative to start)\n",sep=" ")
p <- ggplot(pd,aes(x=timepoint,y=rel_size))
p <- p + geom_vline(xintercept =c(3),color="#fb8072",linetype=1,size=3,alpha=0.3)
p <- p + geom_vline(xintercept =c(2.5,6.5),color="grey30",linetype=2)
p <- p + geom_line(aes(group=ID,color=response),linewidth=0.5,data=pd %>% filter(!is.na(rel_size)))
p <- p + geom_text(aes(label=response,color=response),hjust=0,data=doglabels,nudge_x=0.2,size=2)
p <- p + ggtitle(title)
p <- p + geom_point(aes(color=response,fill=ID,shape=ID),alpha=0.75)
p <- p + scale_y_continuous("lymph node size (relative)")
p <- p + scale_x_continuous("",breaks=timepoints$timepoint,label=timepoints$label,limits=c(0,10),expand = expansion(mult = c(0.05,0.1)))
p <- p + scale_color_manual("dog",values=colors3)
p <- p + scale_fill_manual("dog",values=fills4)
p <- p + scale_shape_manual("dog",values=shapes4)
p <- p + theme_short()
p_lymph_rel <- p

grid <- plot_grid(p_lymph_rel,p_rel_conc,p_rel_TF,p_rel_frag,labels=LETTERS,nrow=2,label_size = 12)
ggsave(grid,filename=pdfname,width=6.5,height=5)