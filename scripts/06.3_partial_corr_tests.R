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
#partial correlation coefficient tests to determine the strength of 
#the relationship between gut microbiota similarity and host phylogenetic
#relatedness, while controlling for host dietary similarity ( %C4 grass in diet)
#and vice-versa
#similar to partial Mantel tests

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load matrix of host divergence times and matrix of
#                   microbiota distances, and make matrix of dietary %C4
################################################################################

load("data/06_divergence_times_matrix.Rdata");
load("data/06_distances_mantel.Rdata");

# make a distance matrix based on host dietary %C4 values
# %C4 reflect the amount of C4 grasses consumed relative to trees,shrubs,leaves
# see manuscript Table S1 for primary sources of these values
diet=data.frame(C4=c(87,2,55,97,41,92,11,92,5,68,91));
rownames(diet)=rownames(divergence_times);
C4.dist=vegdist(diet, method="bray");

################################################################################
#             2. Run partial correlation tests --ALL HERBIVORES
################################################################################
#melt data frames of interest
phylo=melt_dist(divergence_times); phylo$spcomp=paste(phylo$iso1,phylo$iso2, sep="");
diet=melt_dist(as.matrix(C4.dist)); diet$spcomp=paste(diet$iso1,diet$iso2, sep="");
phylo$spcomp==diet$spcomp;

#use for loop to partial correlation coefficient tests on the 4 types of distance matrices
#Bray-Curtis (bray), Jaccard (jac), Weighted Unifrac (wuni), Unweighted Unifrac (unwuni)
mydist=list(mantelbrayall, manteljacall, mantelunifracall, mantelunwunifracall);
names=c("Bray-Curtis", "Jaccard index", "Weighted Unifrac", "Unweighted Unifrac");

for(i in 1:4)
{
  cat(paste("partial correlation coefficient tests, 
              all herbivores, using:", names[i]));
  biota=melt_dist(as.matrix(mydist[[i]])); 
  biota$spcomp=paste(biota$iso1,biota$iso2, sep="");
  #print(biota$spcomp==biota$spcomp);
  final=data.frame(biota$dist, phylo$dist, diet$dist);
  cat("\n examine estimate and p.value for biota-phylo correlation, 
      and biota-diet correlation \n");
  print(pcor(final,method="spearman"));
};


################################################################################
#             3. Run partial correlation tests --BOVIDS ONLY
################################################################################
#subset divergence time matrix to only include Bovids
divergence_times_bov=divergence_times[c(1:4,7:8,10),c(1:4,7:8,10)];

#subset %C4 distance matrix to only include Bovids
C4.dist=as.matrix(C4.dist);
C4.dist.bov=C4.dist[c(1:4,7:8,10),c(1:4,7:8,10)];

#melt divergence times df and diet df
phylo.bov=melt_dist(divergence_times_bov); 
phylo.bov$spcomp=paste(phylo.bov$iso1,phylo.bov$iso2, sep="");

diet.bov=melt_dist(as.matrix(C4.dist.bov)); 
diet.bov$spcomp=paste(diet.bov$iso1,diet.bov$iso2, sep="");

phylo.bov$spcomp==diet.bov$spcomp;

#use for loop to partial correlation coefficient tests on the 4 types of distance matrices
#Bray-Curtis (bray), Jaccard (jac), Weighted Unifrac (wuni), Unweighted Unifrac (unwuni)
mydistbov=list(mantelbraybov, manteljacbov, mantelunifracbov, mantelunwunifracbov);

for(i in 1:4)
{
  cat(paste("partial correlation coefficient tests, 
              Bovids only, using:", names[i]));
  biota2=melt_dist(as.matrix(mydistbov[[i]])); 
  biota2$spcomp=paste(biota2$iso1,biota2$iso2, sep="");
  #print(biota$spcomp==biota$spcomp);
  final2=data.frame(biota2$dist, phylo.bov$dist, diet.bov$dist);
  cat("\n examine estimate and p.value for biota-phylo correlation, 
        and biota-diet correlation \n");
  print(pcor(final2,method="spearman"));
};

