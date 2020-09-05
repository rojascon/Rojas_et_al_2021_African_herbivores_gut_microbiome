#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#MASAI MARA && LAIKIPIA DATA COMBINED
#Code for calculating distance matrices for beta-diversity analyses
#Bray-Curtis, Jaccard, Weighted Unifrac, Unweighted Unifrac
#entire script will run in < 4 minutes

source(file="scripts/background.R"); #load necessary packages and specifications

#load ASV table output by DADA2 && metadata
#samples that did not amplify & baboon samples were removed
#see script "laikipia_mara_get_ASV_table.R"
load("data/laikipia_mara_filtered_ASV_table.Rdata");
load("data/laikipia_mara_filtered_metadata.Rdata");
nrow(metadf)==nrow(asv.tbl);

######################### CALCULATE DISTANCE MATRICES  ##########################
###BRAY-CURTIS distances
bray<-apply(asv.tbl, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distances
jac=(asv.tbl>0)*1;
print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");

###WEIGHTED UNIFRAC distances
#load phylogenetic tree of ASVs needed to calculate Unifrac distances
#see script "laikipia_mara_get_ASV_phylotree_beta.R" for how to generate
load("data/ASV_phylotree_betadiv.Rdata");
wuni<-phyloseq(otu_table(asv.tbl, taxa_are_rows=FALSE), 
               phy_tree(fitGTR$tree));
wuni.dist=UniFrac(wuni, 
                  weighted=TRUE, 
                  normalized=TRUE);

###UNWEIGHTED UNIFRAC distances
unwuni<-phyloseq(otu_table(asv.tbl, taxa_are_rows=FALSE),
                 phy_tree(fitGTR$tree));
unwuni.dist=UniFrac(unwuni, 
                    weighted=FALSE, 
                    normalized=TRUE);

######################### SAVE your distance matrices ##########################
save("bray.dist","jac.dist","wuni.dist","unwuni.dist", 
     file = "data/laikipia_mara_betadiv_dist_objects.Rdata");

