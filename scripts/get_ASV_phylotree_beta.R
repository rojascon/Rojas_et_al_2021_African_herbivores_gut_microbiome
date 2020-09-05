#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for generating phylogenetic tree of ASVs to eventually calculate UNIFRAC distances
#SAMPLES WERE NOT SUBSAMPLED, data as it was output by DADA2 

source(file="scripts/background.R") #load necessary packages and specifications

#load ASV table output by DADA2 R tutorial
#post-removal of samples that did not amplify well (see script "get_ASV_table.R")
load("data/ASV_table.Rdata")

#build phylogenetic tree to calculate Faith's PD 
#these steps are computationally expensive and needs high performance computing
seqs <- getSequences(asv.tbl); 
names(seqs) <- seqs;
alignment <- AlignSeqs(DNAStringSet(seqs), 
                       anchor=NA,verbose=FALSE);
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA"); 
dm <- dist.ml(phangAlign); 
treeNJ <- NJ(dm); 
fit = pml(treeNJ, data=phangAlign);
fitGTR <- update(fit, k=4, inv=0.2); 
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, 
                    optGamma=TRUE, rearrangement = "stochastic", 
                    control = pml.control(trace = 0));

#save the output as .Rdata
save(fitGTR, file="data/ASV_phylotree_betadiv.Rdata");
