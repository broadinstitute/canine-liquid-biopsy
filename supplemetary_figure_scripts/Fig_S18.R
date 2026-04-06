library(tidyverse)
library(cowplot)
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

pdfname <- paste(pdfdir,"Fig_S18.pdf",sep="")
rename_stats <- tibble(order=c(2,1,3),statistic=c("Tumor fraction","cfDNA concentration (ng/mL)","Fragment size ratio"),stat=c("TF","log_DNA_conc","Fragment_size_ratio") ) %>% arrange(order) %>% select(-order)

raw <- as_tibble(read.csv(infile,header=T,sep="\t",na.strings = c("NA","N/A","","UNK")))
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
d <- d %>% mutate(ID=fct_relevel(ID,c("ARM_409177","ARM_412703","ARM_S226485","ARM_408932")))

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
colors4 <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3")


### Make Concentration
pd <- pdAll %>% filter(stat=="log_DNA_conc") 
doglabels <- pd %>% group_by(stat,ID) %>% summarize(timepoint=max(timepoint)) %>% inner_join(pd)
title <- paste("cfDNA concentration\n",sep=" ")

p <- ggplot(pd,aes(x=timepoint,y=value))
p <- p + geom_vline(xintercept =c(3),color="#fb8072",linetype=1,size=3,alpha=0.3)

p <- p + geom_vline(xintercept =c(6.5),color="grey30",linetype=2)
p <- p + geom_line(aes(group=ID,color=ID),linewidth=0.5)
p <- p + geom_text(aes(label=ID,color=ID),hjust=0,data=doglabels,nudge_x=0.2,size=2)
p <- p + ggtitle(title)
p <- p + geom_point(aes(color=ID,shape=ID,y=value),alpha=0.75)
p <- p + scale_y_continuous(str_replace(unique(pd$statistic),"concentration","conc."),breaks=c(-1:3),labels=10**c(-1:3))#expand = expansion(mult = c(0.05, .1)))
p <- p + scale_x_continuous("",breaks=timepoints$timepoint,label=timepoints$label,,expand = expansion(mult = c(0.05,0.1)))
p <- p + scale_color_manual("dog",values=colors4)
p <- p + scale_fill_manual("dog",values=colors4)
p <- p + scale_shape_manual("dog",values=c(15:18))
p <- p + theme_short()
p_conc <- p

### Make TF
pd <- pdAll %>% filter(stat=="TF") 
doglabels <- pd %>% group_by(stat,ID) %>% summarize(timepoint=max(timepoint)) %>% inner_join(pd)
title <- paste("Tumor fraction\n",sep=" ")

p <- ggplot(pd,aes(x=timepoint,y=value))
p <- p + geom_vline(xintercept =c(3),color="#fb8072",linetype=1,size=3,alpha=0.3)

p <- p + geom_vline(xintercept =c(6.5),color="grey30",linetype=2)
p <- p + geom_line(aes(group=ID,color=ID),linewidth=0.5)
p <- p + geom_text(aes(label=ID,color=ID),hjust=0,data=doglabels,nudge_x=0.2,size=2)
p <- p + ggtitle(title)
p <- p + geom_point(aes(color=ID,shape=ID,y=value),alpha=0.75)
p <- p + scale_y_continuous(unique(pd$statistic),limits=c(0,1),breaks=c(0:4)/4) #expand = expansion(mult = c(0.05, .1)))
p <- p + scale_x_continuous("",breaks=timepoints$timepoint,label=timepoints$label,,expand = expansion(mult = c(0.05,0.1)))
p <- p + scale_color_manual("dog",values=colors4)
p <- p + scale_fill_manual("dog",values=colors4)
p <- p + scale_shape_manual("dog",values=c(15:18))
p <- p + theme_short()
p_TF <- p

### Make Fragment size ratio
pd <- pdAll %>% filter(stat=="Fragment_size_ratio") 
doglabels <- pd %>% group_by(stat,ID) %>% summarize(timepoint=max(timepoint)) %>% inner_join(pd)
title <- paste("Fragment size ratio\n",sep=" ")

p <- ggplot(pd,aes(x=timepoint,y=value))
p <- p + geom_vline(xintercept =c(3),color="#fb8072",linetype=1,size=3,alpha=0.3)

p <- p + geom_vline(xintercept =c(6.5),color="grey30",linetype=2)
p <- p + geom_line(aes(group=ID,color=ID),linewidth=0.5)
p <- p + geom_text(aes(label=ID,color=ID),hjust=0,data=doglabels,nudge_x=0.2,size=2)
p <- p + ggtitle(title)
p <- p + geom_point(aes(color=ID,shape=ID,y=value),alpha=0.75)
p <- p + scale_y_continuous(unique(pd$statistic)) #limits=c(0,1),breaks=c(0:4)/4,expand = expansion(mult = c(0.05, .1)))
p <- p + scale_x_continuous("",breaks=timepoints$timepoint,label=timepoints$label,,expand = expansion(mult = c(0.05,0.1)))
p <- p + scale_color_manual("dog",values=colors4)
p <- p + scale_fill_manual("dog",values=colors4)
p <- p + scale_shape_manual("dog",values=c(15:18))
p <- p + theme_short()
p_frag <- p


##### Change relative timepoint 
relative_timepoint <- 0
timepoint_description <- "first timepoint"

relative_timepoint <- 3
timepoint_description <- "treatment timepoint"
pdfname <- str_replace(pdfname,".pdf","_treatment.pdf")

pdRel <- pdAll %>% ungroup() %>% filter(timepoint==relative_timepoint) %>% select(ID,stat,value) %>% rename(firstvalue=value)

pdRelAll <- pdAll %>% left_join(pdRel) %>% mutate(value=value/firstvalue) %>% select(-firstvalue)

### Make Relative Concentration
pd <- pdRelAll %>% filter(stat=="log_DNA_conc") 
doglabels <- pd %>% group_by(stat,ID) %>% summarize(timepoint=max(timepoint)) %>% inner_join(pd)
title <- paste("cfDNA concentration change\n(relative to ",timepoint_description,")",sep="")

p <- ggplot(pd,aes(x=timepoint,y=value))
p <- p + geom_vline(xintercept =c(3),color="#fb8072",linetype=1,size=3,alpha=0.3)

p <- p + geom_vline(xintercept =c(6.5),color="grey30",linetype=2)
p <- p + geom_hline(yintercept=1,color="#535353",linewidth=0.25)

p <- p + geom_line(aes(group=ID,color=ID),linewidth=0.5)
p <- p + geom_text(aes(label=ID,color=ID),hjust=0,data=doglabels,nudge_x=0.2,size=2)
p <- p + ggtitle(title)
p <- p + geom_point(aes(color=ID,shape=ID,y=value),alpha=0.75)
p <- p + scale_y_continuous(str_replace(unique(pd$statistic),"concentration","conc."))
p <- p + scale_x_continuous("",breaks=timepoints$timepoint,label=timepoints$label,,expand = expansion(mult = c(0.05,0.1)))
p <- p + scale_color_manual("dog",values=colors4)
p <- p + scale_fill_manual("dog",values=colors4)
p <- p + scale_shape_manual("dog",values=c(15:18))
p <- p + theme_short()
p_rel_conc <- p


### Make Relative TF
pd <- pdRelAll %>% filter(stat=="TF") 
doglabels <- pd %>% group_by(stat,ID) %>% summarize(timepoint=max(timepoint)) %>% inner_join(pd)
title <- paste("Tumor fraction change\n(relative to ",timepoint_description,")",sep="")

p <- ggplot(pd,aes(x=timepoint,y=value))
p <- p + geom_vline(xintercept=c(3),color="#fb8072",linetype=1,size=3,alpha=0.3)

p <- p + geom_vline(xintercept =c(6.5),color="grey30",linetype=2)
p <- p + geom_hline(yintercept=1,color="#535353",linewidth=0.25)

p <- p + geom_line(aes(group=ID,color=ID),linewidth=0.5)
p <- p + geom_text(aes(label=ID,color=ID),hjust=0,data=doglabels,nudge_x=0.2,size=2)
p <- p + ggtitle(title)
p <- p + geom_point(aes(color=ID,shape=ID,y=value),alpha=0.75)
p <- p + scale_y_continuous(unique(pd$statistic))
p <- p + scale_x_continuous("",breaks=timepoints$timepoint,label=timepoints$label,expand = expansion(mult = c(0.05,0.1)))
p <- p + scale_color_manual("dog",values=colors4)
p <- p + scale_fill_manual("dog",values=colors4)
p <- p + scale_shape_manual("dog",values=c(15:18))
p <- p + theme_short()
p_rel_TF <- p

### Make Relative Fragment size ratio
pd <- pdRelAll %>% filter(stat=="Fragment_size_ratio") 
doglabels <- pd %>% group_by(stat,ID) %>% summarize(timepoint=max(timepoint)) %>% inner_join(pd)
title <- paste("Fragment size ratio change\n(relative to ",timepoint_description,")",sep="")

p <- ggplot(pd,aes(x=timepoint,y=value))
p <- p + geom_vline(xintercept =c(3),color="#fb8072",linetype=1,size=3,alpha=0.3)
p <- p + geom_vline(xintercept =c(6.5),color="grey30",linetype=2)

p <- p + geom_hline(yintercept=1,color="#535353",linewidth=0.25)
p <- p + geom_line(aes(group=ID,color=ID),linewidth=0.5)
p <- p + geom_text(aes(label=ID,color=ID),hjust=0,data=doglabels,nudge_x=0.2,size=2)
p <- p + ggtitle(title)
p <- p + geom_point(aes(color=ID,shape=ID,y=value),alpha=0.75)
p <- p + scale_y_continuous(unique(pd$statistic))
p <- p + scale_x_continuous("",breaks=timepoints$timepoint,label=timepoints$label,,expand = expansion(mult = c(0.05,0.1)))
p <- p + scale_color_manual("dog",values=colors4)
p <- p + scale_fill_manual("dog",values=colors4)
p <- p + scale_shape_manual("dog",values=c(15:18))
p <- p + theme_short()
p_rel_frag <- p

### Make Lymphnodes
pd <- d  %>% select(ID,timepoint,Tumor_size_cm,replicate) %>% distinct() %>% filter(!is.na(Tumor_size_cm)) %>% filter(replicate=="Rep1")
pd <- pd %>% filter(timepoint==0) %>% rename(rel_size=Tumor_size_cm) %>% select(-timepoint) %>% full_join(pd) %>% mutate(rel_size=Tumor_size_cm-rel_size)
doglabels <- pd %>% filter(!is.na(Tumor_size_cm)) %>% group_by(ID) %>% summarize(timepoint=max(timepoint)) %>% inner_join(pd)

title <- paste("Lymph node size\n",sep=" ")
p <- ggplot(pd,aes(x=timepoint,y=Tumor_size_cm))
p <- p + geom_vline(xintercept =c(3),color="#fb8072",linetype=1,size=3,alpha=0.3)

p <- p + geom_vline(xintercept =c(6.5),color="grey30",linetype=2)
p <- p + geom_line(aes(group=ID,color=ID),linewidth=0.5,data=pd %>% filter(!is.na(Tumor_size_cm)))
p <- p + geom_text(aes(label=ID,color=ID,x=timepoint),hjust=0,data=doglabels,nudge_x=0.2,size=2)
p <- p + ggtitle(title)
p <- p + geom_point(aes(color=ID,shape=ID),alpha=0.75)
p <- p + scale_y_continuous("lymph node size")
p <- p + scale_x_continuous("",breaks=timepoints$timepoint,label=timepoints$label,expand = expansion(mult = c(0.05,0.1)))
p <- p + scale_color_manual("dog",values=colors4)
p <- p + scale_fill_manual("dog",values=colors4)
p <- p + scale_shape_manual("dog",values=c(15:18))
p <- p + theme_short()
p_lymph <- p

title <- paste("Lymph node size (relative to start)\n",sep=" ")
p <- ggplot(pd,aes(x=timepoint,y=rel_size))
p <- p + geom_vline(xintercept =c(3),color="#fb8072",linetype=1,size=3,alpha=0.3)

p <- p + geom_vline(xintercept =c(6.5),color="grey30",linetype=2)
p <- p + geom_line(aes(group=ID,color=ID),linewidth=0.5,data=pd %>% filter(!is.na(rel_size)))
p <- p + geom_text(aes(label=ID,color=ID),hjust=0,data=doglabels,nudge_x=0.2,size=2)
p <- p + ggtitle(title)
p <- p + geom_point(aes(color=ID,shape=ID),alpha=0.75)
p <- p + scale_y_continuous("lymph node size (relative)")
p <- p + scale_x_continuous("",breaks=timepoints$timepoint,label=timepoints$label,limits=c(0,10),expand = expansion(mult = c(0.05,0.1)))
p <- p + scale_color_manual("dog",values=colors4)
p <- p + scale_fill_manual("dog",values=colors4)
p <- p + scale_shape_manual("dog",values=c(15:18))
p <- p + theme_short()
p_lymph_rel <- p

### Make x is size lymph node plot
pdAll <- d %>% filter(!is.na(Tumor_size_cm)) %>% select(ID,timepoint,stat,value,replicate,Tumor_size_cm)
pdAll <- pdAll %>% group_by(ID,timepoint,stat,Tumor_size_cm) %>% summarize(value=mean(value))
pdAll <- pdAll %>% left_join(rename_stats)
pdAll <- pdAll %>% mutate(statistic=str_remove(statistic," \\(ng\\/mL\\)"))
pdAll <- pdAll %>% mutate(statistic=str_replace(statistic,"concentration","conc"))
pdAll <- pdAll %>% mutate(statistic=str_remove(statistic," ratio"))
pdAll <- pdAll %>% mutate(statistic=str_replace_all(statistic," ","\n"))
pdAll <- pdAll %>% left_join(timepoints %>% select(timepoint,day))

pdAll$stat <- factor(pdAll$stat,levels=c("log_DNA_conc","TF","Fragment_size_ratio"))
pdAll$statistic <- factor(pdAll$statistic,levels=c("cfDNA\nconc","Tumor\nfraction","Fragment\nsize"))

p <- ggplot(pdAll,aes(x=Tumor_size_cm,y=value))
p <- p + geom_line(aes(group=ID,color=ID),linewidth=1)
p <- p + geom_point(shape=16,size=5,aes(color=ID))
p <- p + geom_point(shape=16,size=4,color="white")

p <- p + geom_text(aes(label=paste("D",day,sep="")),hjust=0.5,vjust=0.5,size=1.75)
p <- p + facet_grid(statistic~ID,scales="free")
p <- p + scale_y_continuous(expand = expansion(mult = c(0.2,0.2)))
p <- p + scale_x_continuous("total lymph node size",expand = expansion(mult = c(0.2,0.2)))
p <- p + scale_color_manual("dog",values=colors4)
p <- p + scale_fill_manual("dog",values=colors4)
p <- p + scale_shape_manual("dog",values=c(15:18))
p <- p + theme_cowplot(12)
p <- p + theme(plot.title = element_blank(),
               axis.title.y=element_blank(),
               axis.title.x = element_text(size=9),axis.text.x = element_text(size=8),
               axis.text.y = element_text(size=7),
               legend.position="none",
               strip.text = element_text(size=6,hjust=0.5,face="bold",color="#f7f7f7"),
               strip.background=element_rect(fill="#525252"),
               panel.border = element_rect(color = "#525252", fill = NA, size = 0.5)
)
p_facet <- p

row2 <- plot_grid(p_facet,labels=LETTERS[9],nrow=1,label_size = 12)
row1 <- plot_grid(p_lymph, p_lymph_rel,labels=LETTERS[1:2],nrow=1,label_size = 12)
grid1 <- plot_grid(p_conc,p_rel_conc,p_TF,p_rel_TF,p_frag,p_rel_frag,ncol=2,label_size = 12,labels=LETTERS[3:8])
grid <- plot_grid(row1,grid1,row2,ncol=1,rel_heights = c(1,3,2))
ggsave(grid,filename=pdfname,width=6.5,height=9)

mainpdfname <- str_replace(pdfname,".pdf",".maintext.pdf")
main_text_grid <- plot_grid(p_lymph_rel,p_rel_conc,p_rel_TF,p_rel_frag,ncol=2,label_size = 12,labels=LETTERS)
ggsave(main_text_grid,filename=mainpdfname,width=6.5,height=4)