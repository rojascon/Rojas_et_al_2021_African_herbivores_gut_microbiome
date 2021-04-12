#################################################################################
#
#                         African herbivores gut microbiome
#                      
#              Rojas et al 2021.Host phylogeny and host ecology structure
#               the mammalian gut microbiota at different taxonomic scales
#
#                               By: Connie Rojas
#                               Created: 5 Sept 2020
#                            Last updated: 18 Jan 2021
#
################################################################################

##CODE FOR: configuring R workspace and printing R version and package versions
#for reader

################################################################################
#             1.  Configure the workspace for subsequent R project scripts                 
################################################################################

#set conditions for R session
rm(list=ls());
options(scipen=999);
options(stringsAsFactors = FALSE) ;

#load necessary packages
library(pacman);
pacman::p_load("car","MASS","dplyr","tidyr","reshape2","vegan","ggplot2",
               "picante","indicspecies","lme4","lmtest","multcomp","grid","phyloseq",
               "Biostrings","QsRutils","phangorn","ape","pheatmap","stringr",
               "dada2","DECIPHER","gridExtra","furrr","harrietr","ggtree",
               "ppcor");

#library(harrietr); library(ggtree); library(ppcor)
################################################################################
#             2. Communicate the R version and package versions to reader 
#                             use SessionInfo()
################################################################################
print("This code was developed with R version 3.6.2");

print("The packages used and their versions were: pheatmap_1.0.12 | phangorn_2.5.5|
QsRutils_0.1.4| Biostrings_2.54.0| phyloseq_1.30.0| multcomp_1.4-15| lmtest_0.9-38| 
lme4_1.1-26| indicspecies_1.7.9| picante_1.8.2| nlme_3.1-151| ape_5.4-1| 
ggplot2_3.3.3| vegan_2.5-7| permute_0.9-5| reshape2_1.4.4| tidyr_1.1.2| dplyr_1.0.3| 
MASS_7.3-53| car_3.0-10| pacman_0.5.1| stringr_1.4.0| DECIPHER_2.14.0| dada2_1.14.1|
gridExtra_2.3| furrr_0.2.1,harrietr_0.2.3,ggtree_2.0.4,ppcor_1.1");



