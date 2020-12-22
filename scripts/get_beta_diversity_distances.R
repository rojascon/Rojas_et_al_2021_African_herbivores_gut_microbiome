#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for calculating distance matrices for beta-diversity analyses
#Bray-Curtis (bray), Jaccard (jac), Weighted Unifrac (wuni), Unweighted Unifrac (unwuni)
#entire script will run in < 4 minutes

source(file="scripts/background.R"); #load necessary packages and specifications
#^^^also loads meta data

#load ASV table output by DADA2 
#samples that did not amplify were removed; see script "get_ASV_table.R"
load("data/ASV_table.Rdata");

#remove baboon samples and samples that did not amplify from metadata file
metadf=meta[meta$Group %in% rownames(asv.tbl),];
metadf=metadf[metadf$species_short!="Baboon",]; 
asvdf=asv.tbl[rownames(asv.tbl) %in% metadf$Group,]; 
nrow(metadf)==nrow(asvdf);

#save formatted metadata file to use in future
metadf=metadf[order(metadf$Group),]; 
metadf$Family=factor(metadf$Family, 
                     levels=c("Bovidae","Giraffidae","Suidae","Equidae","Elephantidae"));
metadf$species_short=factor(metadf$species_short, 
                     levels=c("Buffalo","Cattle","Topi","Eland","Impala","Gazelle",
                              "Dikdik","Giraffe","Zebra","Warthog","Elephant"));
save(metadf, file="data/sample_metadata_beta.Rdata");

######################### ALL HERBIVORES distance matrices  ##########################
###BRAY-CURTIS distance
bray<-apply(asvdf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(asvdf>0)*1;
print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");

###WEIGHTED UNIFRAC distance
#load phylogenetic tree of ASVs needed to calculate Unifrac distances
#see script "get_ASV_phylotree_beta.R" for how to generate
load("data/ASV_phylotree_betadiv.Rdata");
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

######################### BOVIDS ONLY distance matrices  ##########################
#remove non-bovid samples from meta data 
metabovdf=metadf[metadf$Family=="Bovidae",]

#remove non-bovid samples from ASV table
asvbovdf=asvdf[rownames(asvdf) %in% metabovdf$Group,]; 
  
###BRAY-CURTIS distance
braybov<-apply(asvbovdf, 1, function(i) (i/sum(i)));
braybov=as.data.frame(t(braybov));
print(rowSums(braybov));
braybov.dist=vegdist(braybov, method="bray");

###JACCARD distance
jacbov=(asvbovdf>0)*1;
print(rowSums(jacbov));
jacbov.dist=vegdist(jacbov, method="jaccard");

###WEIGHTED UNIFRAC distance
wunibov<-phyloseq(otu_table(asvbovdf, taxa_are_rows=FALSE), 
               phy_tree(fitGTR$tree));
wunibov.dist=UniFrac(wunibov, 
                     weighted=TRUE, 
                     normalized=TRUE);

###UNWEIGHTED UNIFRAC distance
unwunibov<-phyloseq(otu_table(asvbovdf, taxa_are_rows=FALSE),
                 phy_tree(fitGTR$tree));
unwunibov.dist=UniFrac(unwunibov, 
                    weighted=FALSE, 
                    normalized=TRUE);

######################### DIETARY GUILD distance matrices  ##########################
###################these are for making PCoAs for each dietary guild###############
########### one for only grazers, another for only browsers, a third for mixed feeders ####

###GRAZERS
#select only grazer samples from metadata 
metagrazerdf=metadf[metadf$diet_guild=="grazer",];

#select only grazer samples from ASV table
asvgrazerdf=asvdf[rownames(asvdf) %in% metagrazerdf$Group,]; 

###BRAY-CURTIS distance
braygrazer<-apply(asvgrazerdf, 1, function(i) (i/sum(i)));
braygrazer=as.data.frame(t(braygrazer));
print(rowSums(braygrazer));
braygrazer.dist=vegdist(braygrazer, method="bray");

###JACCARD distance
jacgrazer=(asvgrazerdf>0)*1;
print(rowSums(jacgrazer));
jacgrazer.dist=vegdist(jacgrazer, method="jaccard");

###WEIGHTED UNIFRAC distance
wunigrazer<-phyloseq(otu_table(asvgrazerdf, taxa_are_rows=FALSE), 
                  phy_tree(fitGTR$tree));
wunigrazer.dist=UniFrac(wunigrazer, 
                     weighted=TRUE, 
                     normalized=TRUE);

###UNWEIGHTED UNIFRAC distance
unwunigrazer<-phyloseq(otu_table(asvgrazerdf, taxa_are_rows=FALSE),
                    phy_tree(fitGTR$tree));
unwunigrazer.dist=UniFrac(unwunigrazer, 
                       weighted=FALSE, 
                       normalized=TRUE);

###BROWSERS
#select only browser samples from metadata 
metabrowserdf=metadf[metadf$diet_guild=="browser",];

#select only browser samples from ASV table
asvbrowserdf=asvdf[rownames(asvdf) %in% metabrowserdf$Group,]; 

###BRAY-CURTIS distance
braybrowser<-apply(asvbrowserdf, 1, function(i) (i/sum(i)));
braybrowser=as.data.frame(t(braybrowser));
print(rowSums(braybrowser));
braybrowser.dist=vegdist(braybrowser, method="bray");

###JACCARD distance
jacbrowser=(asvbrowserdf>0)*1;
print(rowSums(jacbrowser));
jacbrowser.dist=vegdist(jacbrowser, method="jaccard");

###WEIGHTED UNIFRAC distance
wunibrowser<-phyloseq(otu_table(asvbrowserdf, taxa_are_rows=FALSE), 
                     phy_tree(fitGTR$tree));
wunibrowser.dist=UniFrac(wunibrowser, 
                        weighted=TRUE, 
                        normalized=TRUE);

###UNWEIGHTED UNIFRAC distance
unwunibrowser<-phyloseq(otu_table(asvbrowserdf, taxa_are_rows=FALSE),
                       phy_tree(fitGTR$tree));
unwunibrowser.dist=UniFrac(unwunibrowser, 
                          weighted=FALSE, 
                          normalized=TRUE);

###MIXED-FEEDERS
#select only mixed-feeder samples from meta data 
metamixedfeeddf=metadf[metadf$diet_guild=="mixed_feeder",];

#select only mixed-feedersamples from ASV table
asvmixedfeeddf=asvdf[rownames(asvdf) %in% metamixedfeeddf$Group,]; 

###BRAY-CURTIS distance
braymixedfeed<-apply(asvmixedfeeddf, 1, function(i) (i/sum(i)));
braymixedfeed=as.data.frame(t(braymixedfeed));
print(rowSums(braymixedfeed));
braymixedfeed.dist=vegdist(braymixedfeed, method="bray");

###JACCARD distance
jacmixedfeed=(asvmixedfeeddf>0)*1;
print(rowSums(jacmixedfeed));
jacmixedfeed.dist=vegdist(jacmixedfeed, method="jaccard");

###WEIGHTED UNIFRAC distance
wunimixedfeed<-phyloseq(otu_table(asvmixedfeeddf, taxa_are_rows=FALSE), 
                      phy_tree(fitGTR$tree));
wunimixedfeed.dist=UniFrac(wunimixedfeed, 
                         weighted=TRUE, 
                         normalized=TRUE);

###UNWEIGHTED UNIFRAC distance
unwunimixedfeed<-phyloseq(otu_table(asvmixedfeeddf, taxa_are_rows=FALSE),
                        phy_tree(fitGTR$tree));
unwunimixedfeed.dist=UniFrac(unwunimixedfeed, 
                           weighted=FALSE, 
                           normalized=TRUE);

######################### SAVE your distance matrices ##########################
save("bray.dist", "braybov.dist", "braybrowser.dist", "braygrazer.dist", "braymixedfeed.dist",
     "jac.dist","jacbov.dist", "jacbrowser.dist", "jacgrazer.dist", "jacmixedfeed.dist", 
     "unwuni.dist", "unwunibov.dist", "unwunibrowser.dist", "unwunigrazer.dist", "unwunimixedfeed.dist",
     "wuni.dist", "wunibov.dist", "wunibrowser.dist", "wunigrazer.dist", "wunimixedfeed.dist", 
     file = "data/betadiv_dist_objects.Rdata");
