#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#MASAI MARA && LAIKIPIA DATA COMBINED
#Code for beta-diversity analyses (PERMANOVA)

source(file="scripts/background.R"); #load necessary packages and specifications

#load metadata with baboon samples removed, as well as samples that did not amplify
load("data/laikipia_mara_filtered_metadata.Rdata");

#load distance matrices computed from ASV abundances
#see script "laikipia_mara_get_beta_diversity_distances.R"
load("data/laikipia_mara_betadiv_dist_objects.Rdata");

#######################  PERMANOVAS  ##############################
############## GLOBAL MODEL
print("BRAY/JACCARD: PERMANOVA all host predictors");
print(adonis(bray.dist ~     ##bray.dist, jac.dist
               diet_guild+
               species_short+    
               region+
               sample_month,
             data=metadf,
             method = "bray",     #bray, jaccard
             permutations = 999));

print("UNIFRAC: PERMANOVA all host predictors");
print(adonis(wuni.dist ~     ##wuni.dist, unwuni.dist
               diet_guild+
               species_short+    
               region+
               sample_month,
             data=metadf,
             permutations = 999));

############## WITHIN A SPECIES, do samples vary by region?
print("BRAY/JACCARD: PERMANOVA within host species");
print(adonis(bray.dist ~     ##bray.dist, jac.dist
               region+
               sample_month, 
             strata=metadf$species_short,
             data=metadf,
             method = "bray",    #bray, jaccard
             permutations = 999));

print("UNIFRAC: PERMANOVA within host species");
print(adonis(wuni.dist ~     ##wuni.dist, unwuni.dist
               region+
               sample_month, 
             strata=metadf$species_short,
             data=metadf,
             permutations = 999));

############## WITHIN A DIET GUILD, do samples vary by species and region?
print("BRAY/JACCARD: PERMANOVA within dietary guild");
print(adonis(bray.dist ~     ##bray.dist, jac.dist
               species_short+
               region+
               sample_month,
             strata=metadf$diet_guild,
             data=metadf,
             method = "bray",    #bray, jaccard
             permutations = 999));

print("UNIFRAC: PERMANOVA within dietary guild");
print(adonis(wuni.dist ~     ##wuni.dist, unwuni.dist
               species_short+
               region+
               sample_month,
             strata=metadf$diet_guild,
             data=metadf,
             permutations = 999));
