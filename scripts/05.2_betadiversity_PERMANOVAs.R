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

##CODE FOR: conducting PERMANOVA tests based on 4 distance matrices:
#Bray-Curtis (bray), Jaccard (jac), Weighted Unifrac (wuni), Unweighted Unifrac (unwuni)

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load beta diversity distance matrices (e.g. dissimilarity matrices) 
#                         & filtered metadata table                 
################################################################################
load("data/01_sample_metadata_filtered.Rdata");
load("data/05_dissimilarity_distances_beta.Rdata")


################################################################################
#             2. conduct PERMANOVAs -- ALL HERBIVORES (host family)                
################################################################################
#basically asking, does gut microbiota structure vary with sample month, host
#dietary guild, and host family?
#because the same test is done on each type of distance matrix [4 total], we
#will use a for loop

#set variables for for loop
mydist=list(bray.dist, jac.dist, wuni.dist, unwuni.dist)
names=c("Bray-Curtis", "Jaccard index", "Weighted Unifrac", "Unweighted Unifrac")
met=c("bray","jaccard","bray","jaccard") 

#run for loop to conduct the 4 PERMANOVA tests
for(i in 1:4)
{
  print(paste("PERMANOVA test, all herbivores, using:", names[i]));
  print(adonis(mydist[[i]]~     
                 sample_month+
                 diet_guild+
                 Family,
               data=meta,
               method = met[i],     
               permutations = 999));
};

################################################################################
#             3. conduct PERMANOVAs -- ALL HERBIVORES (host species)                
################################################################################
#basically asking, does gut microbiota structure vary with sample month, host
#dietary guild, and host species?
#because the same test is done on each type of distance matrix [4 total], we
#will use a for loop

#set variables for for loop
mydist=list(bray.dist, jac.dist, wuni.dist, unwuni.dist)
names=c("Bray-Curtis", "Jaccard index", "Weighted Unifrac", "Unweighted Unifrac")
met=c("bray","jaccard","bray","jaccard") 

#run for loop to conduct the 4 PERMANOVA tests
for(i in 1:4)
{
  print(paste("PERMANOVA test, all herbivores, using:", names[i]));
  print(adonis(mydist[[i]]~     
                 sample_month+
                 diet_guild+
                 species_short,   
               data=meta,
               method = met[i],     
               permutations = 999));
};

################################################################################
#             4. conduct PERMANOVAs -- BOVIDS ONLY                
################################################################################
#same tests as above, but restricting dataset to bovids only

#remove non-bovid samples from dataset
metabov=meta[meta$Family=="Bovidae",]
metabov$species_short=factor(metabov$species_short, 
                               levels=c("Buffalo","Cattle","Topi","Eland","Impala",
                                        "Gazelle","Dikdik"));
#set variables for for loop
mydist=list(braybov.dist, jacbov.dist, wunibov.dist, unwunibov.dist)
#names and met variables are same as above

#run for loop to conduct the 4 PERMANOVA tests
for(i in 1:4)
{
  print(paste("PERMANOVA test, Bovids only, using:", names[i]));
  print(adonis(mydist[[i]]~     
                 sample_month+
                 diet_guild+
                 species_short,   
               data=metabov,
               method = met[i],      
               permutations = 999));
};

