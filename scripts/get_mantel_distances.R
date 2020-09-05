#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for calculating average microbiota distances (Bray-Curtis, Jaccard, Weighted Unifrac, 
#Unweighted Unifrac) for MANTEL tests 

#SIMILAR BUT NOT IDENTICAL TO "get_beta_diversity_distances.R"
#HERE WE ARE AVERAGING MICROBIOTA PROFILES FOR EACH HOST HERBIVORE SPECIES SO THAT WE HAVE
#11 x 11 MATRICES INSTEAD OF 165 x 165 MATRICES

source(file="scripts/background.R"); #load necessary packages and specifications

#load metadata with baboon samples removed, as well as samples that did not amplify
load("data/sample_metadata_beta.Rdata");

#load ASV table output by DADA2 
#samples that did not amplify were removed; see script "get_ASV_table.R"
load("data/ASV_table.Rdata");

#remove baboon samples from ASV table
asvdf=asv.tbl[rownames(asv.tbl) %in% metadf$Group,]; 
nrow(metadf)==nrow(asvdf);

##transform ASV counts to ASV proportions
asvman<- apply(asvdf, 1, function(i) (i/sum(i)))
colSums(asvman) 

#replace column names with their host species names'
colnames(asvman)=metadf$species_short[
  match(colnames(asvman),metadf$Group)];

#average the microbiota profiles for each host species
myspecies=unique(metadf$species_short)
col_names=colnames(asvman); 

for(i in 1:length(myspecies))
{
  A=grep(pattern=myspecies[i], x=col_names);
  B=rowMeans(asvman[,A]);
  assign(paste(myspecies[i]), B); 
  }

#make master data frame of the averaged microbiota profiles
profiles=cbind(Cattle, Dikdik, Impala, Topi, Elephant, Zebra, Eland, Buffalo, 
               Giraffe,Gazelle, Warthog)

tprofiles=t(profiles); 
rowSums(tprofiles);

#calculate distance matrices
#do it separately for all herbivores and for bovids only

######################### ALL HERBIVORES mantel distance matrices  ##########################
####BRAY-CURTIS distances
mantelbrayall=vegdist(tprofiles, method="bray"); 

####JACCARD distances
j=(tprofiles>0)*1
manteljacall=vegdist(j, method="jaccard");

###WEIGHTED UNIFRAC distances
#load phylogenetic tree of ASVs needed to calculate Unifrac distances
#see script "get_phylotree_beta.R" for how to generate
load("data/ASV_phylotree_betadiv.Rdata");
uni<-phyloseq(otu_table(tprofiles, taxa_are_rows=FALSE),
              phy_tree(fitGTR$tree));
mantelunifracall=UniFrac(uni, 
                    weighted=TRUE, 
                    normalized=TRUE);

###UNWEIGHTED UNIFRAC distances
mantelunwunifracall=UniFrac(uni, 
                         weighted=FALSE, 
                         normalized=TRUE);

######################### BOVIDS ONLY distance matrices  ##########################
#remove non-bovid samples from averaged ASV table
tbovids=tprofiles[c(1:4,7:8,10),]

####BRAY-CURTIS distances
mantelbraybov=vegdist(tbovids, method="bray"); 

####JACCARD distances
j=(tbovids>0)*1
manteljacbov=vegdist(j, method="jaccard");

###WEIGHTED UNIFRAC distances
#load phylogenetic tree of ASVs needed to calculate Unifrac distances
#see script "get_phylotree_beta.R" for how to generate
load("data/ASV_phylotree_betadiv.Rdata");
uni<-phyloseq(otu_table(tbovids, taxa_are_rows=FALSE),
              phy_tree(fitGTR$tree));
mantelunifracbov=UniFrac(uni, 
                         weighted=TRUE, 
                         normalized=TRUE);

###UNWEIGHTED UNIFRAC distances
mantelunwunifracbov=UniFrac(uni, 
                            weighted=FALSE, 
                            normalized=TRUE);

######################### SAVE your mantel matrices ##########################
save(mantelbrayall, mantelbraybov, manteljacall, manteljacbov, mantelunifracall, 
     mantelunifracbov, mantelunwunifracall, mantelunwunifracbov, 
    file = "data/mantel_dist_objects.Rdata")
