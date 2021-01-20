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

##CODE FOR: running 
#A)mantel tests of host microbiota dissimilarity vs host divergence time 

#B) partial mantel tests of host microbiota dissimilarity vs host divergence time,
#while controlling for dietary similarity 

source(file="scripts/00_background.R"); #load necessary packages and specifications

################################################################################
#             1. Load matrix of host divergence times and matrix of
#                   microbiota distances
################################################################################

load("data/06_divergence_times_matrix.Rdata");
load("data/06_distances_mantel.Rdata");

################################################################################
#             2. Run Mantel test --ALL HERBIVORES
################################################################################

#use for loop to run mantel test on the 4 types of distance matrices
#Bray-Curtis (bray), Jaccard (jac), Weighted Unifrac (wuni), Unweighted Unifrac (unwuni)
mydist=list(mantelbrayall, manteljacall, mantelunifracall, mantelunwunifracall);
names=c("Bray-Curtis", "Jaccard index", "Weighted Unifrac", "Unweighted Unifrac");

for(i in 1:4)
{
  print(paste("Mantel test, all herbivores, using:", names[i]));
  print(mantel(mydist[[i]], divergence_times, 
               method="spear", 
               permutations=999));
};
 
################################################################################
#             3. Run Mantel test --BOVIDS ONLY
################################################################################
#subset divergence time matrix to only bovids
divergence_times_bov=divergence_times[c(1:4,7:8,10),c(1:4,7:8,10)];

#use for loop to run mantel test on each distance matrix
mydistbov=list(mantelbraybov, manteljacbov, mantelunifracbov, mantelunwunifracbov);

for(i in 1:4)
{
  print(paste("Mantel test, all herbivores, using:", names[i]));
  print(mantel(mydistbov[[i]], divergence_times_bov, 
               method="spear", 
               permutations=999));
};
################################################################################
#             4. Run Partial Mantel test --ALL HERBIVORES
################################################################################

# make a distance matrix based on host dietary %C4 values
# %C4 reflect the amount of trees, leaves, and forbs consumed relative to grasses
# see manuscript Table S4 for primary sources
diet=data.frame(C4=c(87,2,55,97,41,92,11,92,5,68,91));
rownames(diet)=rownames(divergence_times);
C4.dist=vegdist(diet, method="bray");

#use for loop to run partial mantel test on each distance matrix
for(i in 1:4)
{
  print(paste("Partial Mantel test, all herbivores, using:", names[i]));
  print(mantel.partial(mydist[[i]], divergence_times,
               C4.dist,
               method="spear", 
               permutations=999));
};


################################################################################
#             4. Run Partial Mantel test --BOVIDS ONLY
################################################################################

#subset %C4 distance matrix to only include Bovids
C4.dist=as.matrix(C4.dist);
C4.dist.bov=C4.dist[c(1:4,7:8,10),c(1:4,7:8,10)];

#use for loop to run partial mantel test on each distance matrix
for(i in 1:4)
{
  print(paste("Partial Mantel test, all herbivores, using:", names[i]));
  print(mantel.partial(mydistbov[[i]], divergence_times_bov,
                       C4.dist.bov,
                       method="spear", 
                       permutations=999));
};
