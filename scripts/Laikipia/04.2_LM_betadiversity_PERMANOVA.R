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

##Note these are analyses for the COMBINED LAIKIPIA AND MASAI MARA DATASET
#working from Laikipia subdirectories within the main directories

##CODE FOR: conducting PERMANOVA tests based on 4 distance matrices:
#Bray-Curtis (bray), Jaccard (jac), Weighted Unifrac (wuni), Unweighted Unifrac (unwuni)

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load beta diversity distance matrices (e.g. dissimilarity matrices) 
#                         & filtered metadata table                 
################################################################################
load("data/Laikipia/01_LM_sample_metadata_filtered.Rdata");
load("data/Laikipia/04_LM_dissimilarity_distances_beta.Rdata")


################################################################################
#             2. conduct PERMANOVAs                  
################################################################################
#basically asking, does gut microbiota structure vary with sample month, host
#dietary guild, host geographic region, and host species?
#because the same test is done on each type of distance matrix [4x total], we
#will use a for loop

#set variables for for loop
mydist=list(bray.dist, jac.dist, wuni.dist, unwuni.dist)
names=c("Bray-Curtis", "Jaccard index", "Weighted Unifrac", "Unweighted Unifrac")
met=c("bray","jaccard","bray","jaccard") 

#run for loop to conduct the 4 PERMANOVA tests
for(i in 1:4)
{
  print(paste("PERMANOVA test, mara & laikipia herbivores, using:", names[i]));
  print(adonis(mydist[[i]]~     
                 sample_month+
                 region+
                 diet_guild+
                 species_short,   
               data=meta,
               method = met[i],     
               permutations = 999));
};

