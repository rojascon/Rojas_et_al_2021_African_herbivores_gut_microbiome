#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#MASAI MARA && LAIKIPIA DATA COMBINED

##Code for testing phylosymbiosis (mantel tests of host divergence time vs. 
#host microbiota dissimilarity) 

source(file="scripts/background.R"); #load necessary packages and specifications

#load metadata with baboon samples removed, as well as samples that did not amplify
load("data/laikipia_mara_filtered_metadata.Rdata");

#load host species divergence times
#see script "get_host_divergence_times.R"
load("data/host_divergence_times.Rdata");

#load gut microbiota distance matrices built for this mantel test
#see script "laikipia_mara_get_mantel_distances.R"
#one set of distance matrices for all herbivores; another set for bovids only
load("data/LM_mantel_dist_objects.Rdata");

##include only host species of interest in table of divergence times
LMdivergence_times=divergence_times[c(1,3,5:9,11),c(1,3,5:9,11)]

#possible distance matries are
#Bray-Curtis all host species 
######################## RUN MANTEL TESTS ##########################

#run mantel test---all herbivores
#possible distance matrices are: BRAY-CURTIS (LMmantelbrayall), JACCARD (LMmanteljacall),
#WEIGHTED-UNIFRAC (LMmantelunifracall), UNWEIGHTED-UNIFRAC (LMmantelunwunifracall)

colnames(LMdivergence_times)==
  colnames(as.matrix(LMmantelbrayall)); #pick distance matrix here

mantel(LMdivergence_times, LMmantelunifracall, #pick distance matrix here
       method="spear", 
       permutations=999); 

#run mantel test---bovids only
#possible distance matrices are: BRAY-CURTIS (LMmantelbraybov), JACCARD (LMmanteljacbov),
#WEIGHTED-UNIFRAC (LMmantelunifracbov), UNWEIGHTED-UNIFRAC (LMmantelunwunifracbov)

#include only bovids in table of divergence times
divergence_times_bov=LMdivergence_times[c(1:2,5:6),c(1:2,5:6)];

colnames(divergence_times_bov)==
  colnames(as.matrix(LMmantelbraybov));   #pick distance matrix here

mantel(divergence_times_bov, LMmantelunwunifracbov,   #pick distance matrix here
       method="spear", 
       permutations=999);
