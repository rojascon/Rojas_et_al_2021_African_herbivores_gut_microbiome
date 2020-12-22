#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#MASAI MARA && LAIKIPIA DATA COMBINED

#Code for calculating microbiota distances (Bray-Curtis, Jaccard, Weighted Unifrac, 
#Unweighted Unifrac) from averaged microbiota profiles for MANTEL tests 

#SIMILAR BUT NOT IDENTICAL TO "get_beta_diversity_distances.R"
#HERE WE ARE AVERAGING MICROBIOTA PROFILES FOR EACH HOST SPECIES SO THAT WE HAVE
#11 x 11 MATRICES INSTEAD OF 165 x 165 MATRICES

source(file="scripts/background.R"); #load necessary packages and specifications

#load ASV table output by DADA2 && metadata
#samples that did not amplify & baboon samples were removed
#see script "laikipia_mara_get_ASV_table.R"
load("data/laikipia_mara_filtered_ASV_table.Rdata");
load("data/laikipia_mara_filtered_metadata.Rdata");
nrow(metadf)==nrow(asv.tbl);

##transform ASV counts to ASV proportions
asvman=as.matrix(t(asv.tbl));
asvman<- apply(asv.tbl, 1, function(i) (i/sum(i)));
colSums(asvman);

#replace column names with their host species names'
colnames(asvman)=metadf$species_short[
  match(colnames(asvman),metadf$Group)];

#average the microbiota profiles for each host species
#so all of the buffalo samples get averaged into 1 buffalo profile
myspecies=unique(metadf$species_short);
col_names=colnames(asvman); 

for(i in 1:length(myspecies))
{
  A=grep(pattern=myspecies[i], x=col_names);
  B=rowMeans(asvman[,A]);
  assign(paste(myspecies[i]), B); 
};

#make master data frame of the averaged microbiota profiles
profiles=cbind(Cattle, Impala, Elephant, Zebra, Eland, Buffalo, 
               Giraffe,Warthog);

tprofiles=t(profiles); 
rowSums(tprofiles);

#calculate distance matrices
#do it separately for all herbivores and for bovids only

######################### ALL HERBIVORES mantel distance matrices  ##########################
####BRAY-CURTIS distance
LMmantelbrayall=vegdist(tprofiles, method="bray"); 

####JACCARD distance
j=(tprofiles>0)*1
LMmanteljacall=vegdist(j, method="jaccard");

###WEIGHTED UNIFRAC distance
#load phylogenetic tree of ASVs needed to calculate Unifrac distances
#see script "get_phylotree_beta.R" for how to generate
load("data/ASV_phylotree_betadiv.Rdata");
uni<-phyloseq(otu_table(tprofiles, taxa_are_rows=FALSE),
              phy_tree(fitGTR$tree));
LMmantelunifracall=UniFrac(uni, 
                         weighted=TRUE, 
                         normalized=TRUE);

###UNWEIGHTED UNIFRAC distance
LMmantelunwunifracall=UniFrac(uni, 
                            weighted=FALSE, 
                            normalized=TRUE);

######################### BOVIDS ONLY distance matrices  ##########################
#remove non-bovid samples from ASV table
tbovids=tprofiles[c(1:2,5:6),]

####BRAY-CURTIS distance
LMmantelbraybov=vegdist(tbovids, method="bray"); 

####JACCARD distance
j=(tbovids>0)*1
LMmanteljacbov=vegdist(j, method="jaccard");

###WEIGHTED UNIFRAC distance
#load phylogenetic tree of ASVs needed to calculate Unifrac distances
#see script "get_phylotree_beta.R" for how to generate
load("data/ASV_phylotree_betadiv.Rdata");
uni<-phyloseq(otu_table(tbovids, taxa_are_rows=FALSE),
              phy_tree(fitGTR$tree));
LMmantelunifracbov=UniFrac(uni, 
                         weighted=TRUE, 
                         normalized=TRUE);

###UNWEIGHTED UNIFRAC distance
LMmantelunwunifracbov=UniFrac(uni, 
                            weighted=FALSE, 
                            normalized=TRUE);

######################### SAVE your mantel matrices ##########################
save(LMmantelbrayall, LMmantelbraybov, LMmanteljacall, LMmanteljacbov, LMmantelunifracall, 
     LMmantelunifracbov, LMmantelunwunifracall, LMmantelunwunifracbov, 
     file = "data/LM_mantel_dist_objects.Rdata");

