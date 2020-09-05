#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for beta-diversity analyses (PERMANOVA)

source(file="scripts/background.R"); #load necessary packages and specifications

#load metadata with baboon samples removed, as well as samples that did not amplify
load("data/sample_metadata_beta.Rdata");

#load distance matrices computed from ASV abundances
#see script "get_beta_diversity_distances.R"
#one set of distance matrices for all herbivores; another set for bovids only; 
#a third set for grazers, browsers, mixed feeders
load("data/betadiv_dist_objects.Rdata")

#######################  PERMANOVAS  ##############################
############## across all herbivores
print("BRAY/JACCARD: PERMANOVA all host predictors, all herbivores");
print(adonis(bray.dist~diet_guild+Family+species_short+   ##bray.dist, jac.dist
         sample_month+rain_two_weeks+tmax_two_weeks+
         tmin_two_weeks,
       data=metadf,
       method = "bray",      #bray, jaccard
       permutations = 999));


print("UNIFRAC: PERMANOVA all host predictors, all herbivores"); 
print(adonis(wuni.dist~diet_guild+Family+species_short+    ##wuni.dist, unwuni.dist
         sample_month+rain_two_weeks+tmax_two_weeks+
         tmin_two_weeks,
       data=metadf,
       permutations = 999));

############## only bovids
#make a bovid only metadata file 
metabovdf=metadf[metadf$Family=="Bovidae",]

print("BRAY/JACCARD: PERMANOVA most host predictors, only bovids");
print(adonis(braybov.dist~diet_guild+species_short+sample_month+  ##braybov.dist, jacbov.dist
         rain_two_weeks+tmax_two_weeks+tmin_two_weeks, 
       data=metabovdf, 
       method = "bray",        #bray, jaccard
       permutations = 999));

print("UNIFRAC: PERMANOVA most host predictors, only bovids");
print(adonis(wunibov.dist~diet_guild+species_short+sample_month+  ##wunibov.dist, unwunibov.dist
         rain_two_weeks+tmax_two_weeks+tmin_two_weeks, 
       data=metabovdf, 
       permutations = 999));

############### all herbivores, controlling for dietary guild
print("BRAY/JACCARD: PERMANOVA controlling for dietary guild, all herbivores");
print(adonis(bray.dist~species_short+sample_month+rain_two_weeks+  ##bray.dist, jac.dist
         tmax_two_weeks+tmin_two_weeks, 
       strata=metadf$diet_guild, 
       data=metadf, 
       method = "bray",         #bray, jaccard
       permutations = 999));

print("UNIFRAC: PERMANOVA controlling for dietary guild, all herbivores");
print(adonis(wuni.dist~species_short+sample_month+rain_two_weeks+  ##wuni.dist, unwuni.dist
         tmax_two_weeks+tmin_two_weeks, 
       strata=metadf$diet_guild, 
       data=metadf,
       permutations = 999));

