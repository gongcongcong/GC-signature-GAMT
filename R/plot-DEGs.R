library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(UpSetR)
# load data ---------------------------------------------------------------

dat <- readxl::read_xlsx('data/DE_RNA_20180405.xlsx',sheet='sum') %>%
        gather(key='group',
               value='gene') 
dat <- dat[complete.cases({dat}$gene),]
table(dat$group)
dat_list <- list('t1_vs_t2'=unlist(dat[dat$group=='t1_vs_t2','gene']),'t1_vs_t3'=unlist(dat[dat$group=='t1_vs_t3','gene']),'t1_vs_t4'=unlist(dat[dat$group=='t1_vs_t4','gene']),'t2_vs_t3'=unlist(dat[dat$group=='t2_vs_t3','gene']),'t2_vs_t4'=unlist(dat[dat$group=='t2_vs_t4','gene']),'t3_vs_t4'=unlist(dat[dat$group=='t3_vs_t4','gene']))
dat_for_upset <- UpSetR::fromList(dat_list)
p_upset_intersects<- upset(dat_for_upset,order.by='freq')
pdf(file='plots/p_upset_DEGs_intersects.pdf')
p_upset_intersects
dev.off()
