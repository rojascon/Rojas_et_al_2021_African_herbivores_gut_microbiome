#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#MASAI MARA && LAIKIPIA DATA COMBINED
#Code for beta-diversity analyses (PERMANOVAs)

source(file="scripts/background.R"); #load necessary packages and specifications

#load metadata with baboon samples removed, as well as samples that did not amplify
load("data/laikipia_mara_filtered_metadata.Rdata");

#load distance matrices computed from ASV abundances
#see script "laikipia_mara_get_beta_diversity_distances.R"
load("data/laikipia_mara_betadiv_dist_objects.Rdata");

#######################  RUN PERMANOVAS  ##############################
print("BRAY/JACCARD: PERMANOVA all host predictors");
print(adonis(bray.dist ~     ##bray.dist, jac.dist
               sample_month+
               region+
               diet_guild+
               species_short,
             data=metadf,
             method = "bray",     #bray, jaccard
             permutations = 999));

print("UNIFRAC: PERMANOVA all host predictors");
print(adonis(unwuni.dist ~     ##wuni.dist, unwuni.dist
               sample_month+
               region+
               diet_guild+
               species_short,
             data=metadf,
             permutations = 999));

