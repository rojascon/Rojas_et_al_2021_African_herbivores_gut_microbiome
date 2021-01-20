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

##CODE FOR: calculating distance matrices for beta-diversity analyses
#the 4 types of distance matrices are:
#Bray-Curtis (bray), Jaccard (jac), Weighted Unifrac (wuni), Unweighted Unifrac (unwuni)

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load filtered ASV abundance table, ASV taxonomy table and 
#                         filtered metadata table                 
################################################################################
tax=read.csv("data/Laikipia/00_LM_ASV_taxonomy_silva_v132.csv", header=T, row.names=1);
load("data/Laikipia/01_LM_ASV_table_filtered.Rdata");
load("data/Laikipia/01_LM_sample_metadata_filtered.Rdata");


################################################################################
#             2. Create phylogenetic tree of ASV sequences
#     this is necessary for calculating Unifrac distances later
################################################################################
##these steps are computationally expensive and needs high performance computing
#that is why they are commented out
#transfer this part of the code and the necessary files to R on the HPC at your University

# seqs <- getSequences(as.matrix(asvdf)); 
# names(seqs) <- seqs;
# alignment <- AlignSeqs(DNAStringSet(seqs), 
#                        anchor=NA,verbose=FALSE);
# phangAlign <- phyDat(as(alignment, "matrix"), type="DNA"); 
# dm <- dist.ml(phangAlign); 
# treeNJ <- NJ(dm); 
# fit = pml(treeNJ, data=phangAlign);
# fitGTR <- update(fit, k=4, inv=0.2); 
# fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, 
#                     optGamma=TRUE, rearrangement = "stochastic", 
#                     control = pml.control(trace = 0));

#save the output as an .Rdata file
#save(fitGTR, file="data/Laikipia/04_LM_ASV_phylotree.Rdata");


################################################################################
#             3. Generate the 4 types of distance matrices
################################################################################
###BRAY-CURTIS distance
bray<-apply(asvdf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
#print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(asvdf>0)*1;
#print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");

###WEIGHTED UNIFRAC distance
load("data/Laikipia/04_LM_ASV_phylotree.Rdata");
wuni<-phyloseq(otu_table(asvdf, taxa_are_rows=FALSE), 
               phy_tree(fitGTR$tree));
wuni.dist=UniFrac(wuni, 
                  weighted=TRUE, 
                  normalized=TRUE);

###UNWEIGHTED UNIFRAC distance
unwuni<-phyloseq(otu_table(asvdf, taxa_are_rows=FALSE),
                 phy_tree(fitGTR$tree));
unwuni.dist=UniFrac(unwuni, 
                    weighted=FALSE, 
                    normalized=TRUE);


################################################################################
#             4. Save all distance matrices as a single R.data object
################################################################################

save("bray.dist","jac.dist","wuni.dist","unwuni.dist", 
     file = "data/laikipia/04_LM_dissimilarity_distances_beta.Rdata");