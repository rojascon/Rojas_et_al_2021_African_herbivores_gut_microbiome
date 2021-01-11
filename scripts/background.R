#set conditions for R session
rm(list=ls());
options(scipen=999);
library(pacman);
pacman::p_load("car","MASS","plyr","dplyr","tidyr","reshape2","data.table","stringr",
               "vegan","ggplot2","RColorBrewer","gplots","PMCMR","MuMIn","dada2",
               "DECIPHER","picante","indicspecies","nlme","lme4","lmerTest","lmtest",
               "multcomp","emmeans","e1071","glmmTMB","betareg","gdata","grid",
               "gridExtra","lattice","phyloseq","Biostrings","ranacapa","QsRutils",
               "phangorn", "DESeq2","ape","furrr","future","pheatmap");


#format the order of levels for factors of host predictors
meta=read.csv("data/sample_metadata.csv", header=T);

meta$sample_month=factor(meta$sample_month, 
                         levels=c("march","april","may","june"));
meta$Order=factor(meta$Order, 
                  levels=c("Artiodactyla","Perissodactyla","Proboscidea","Primata"));
meta$Family=factor(meta$Family, 
                   levels=c("Bovidae","Giraffidae","Suidae","Equidae","Elephantidae","Cercopithecidae"));
meta$species_short=factor(meta$species_short, 
                          levels=c("Buffalo","Cattle","Topi","Eland","Impala","Gazelle",
                                                       "Dikdik","Giraffe","Zebra","Warthog","Elephant","Baboon"))
meta$diet_guild=factor(meta$diet_guild, 
                       levels=c("grazer","browser","mixed_feeder"));

#check version of your loaded packages
#packageVersion('ape')
