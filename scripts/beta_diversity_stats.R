#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for beta-diversity analyses (PERMANOVAs)

source(file="scripts/background.R"); #load necessary packages and specifications

#load metadata with baboon samples removed, as well as samples that did not amplify
load("data/sample_metadata_beta.Rdata");

#load distance matrices computed from ASV abundances
#see script "get_beta_diversity_distances.R"
#one set of distance matrices for all herbivores; another set for bovids only; 
#a third set for only grazers, only browsers, only mixed feeders
load("data/betadiv_dist_objects.Rdata")

#######################  PERMANOVAS  ##############################
############## across all herbivores
print("BRAY/JACCARD: PERMANOVA all host predictors, all herbivores");
print(adonis(bray.dist~     ##bray.dist, jac.dist
               sample_month+
               diet_guild+
               Family+
               species_short,   
                data=metadf,
                method = "bray",      #bray, jaccard
                permutations = 999));


print("UNIFRAC: PERMANOVA all host predictors, all herbivores"); 
print(adonis(wuni.dist~              ##wuni.dist, unwuni.dist
               sample_month+
               diet_guild+
               Family+
               species_short,    
                data=metadf,
                permutations = 999));


############## only bovids
#make a bovid only metadata file 
metabovdf=metadf[metadf$Family=="Bovidae",]
metabovdf$species_short=factor(metabovdf$species_short, 
                            levels=c("Buffalo","Cattle","Topi","Eland","Impala",
                                     "Gazelle","Dikdik"));

print("BRAY/JACCARD: PERMANOVA most host predictors, only bovids");
print(adonis(jacbov.dist~               ##braybov.dist, jacbov.dist
               sample_month+
               diet_guild+
               species_short,  
                data=metabovdf, 
                method = "jaccard",        #bray, jaccard
                permutations = 999));

print("UNIFRAC: PERMANOVA most host predictors, only bovids");
print(adonis(unwunibov.dist~           ##wunibov.dist, unwunibov.dist
               sample_month+
               diet_guild+
               species_short,  
                data=metabovdf, 
                permutations = 999));

############### all herbivores, controlling for dietary guild and local environmental conditions
print("BRAY/JACCARD: PERMANOVA controlling for dietary guild, all herbivores");
print(adonis(jac.dist~                 ##bray.dist, jac.dist
               sample_month+
               species_short,  
             strata=metadf$diet_guild, 
             data=metadf, 
             method = "jaccard",         #bray, jaccard
             permutations = 999));

print("UNIFRAC: PERMANOVA controlling for dietary guild, all herbivores");
print(adonis(unwuni.dist~               ##wuni.dist, unwuni.dist
              sample_month+
              species_short,  
            strata=metadf$diet_guild, 
            data=metadf,        
            permutations = 999));


