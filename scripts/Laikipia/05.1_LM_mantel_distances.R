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

##CODE FOR: generating microbiota distance matrices Bray-Curtis (bray), 
#Jaccard (jac), Weighted Unifrac (wuni), Unweighted Unifrac (unwuni) for Mantel tests

#These distances are based on microbiota profiles that were AVERAGED across host species 
#Their dimensions match the dimensions of the matrix of host divergence times

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load filtered ASV table & filtered metadata table                 
################################################################################
load("data/Laikipia/01_LM_sample_metadata_filtered.Rdata");
load("data/Laikipia/01_LM_ASV_table_filtered.Rdata");

#turn ASV counts to ASV proportions
asvman<- apply(asvdf, 1, function(i) (i/sum(i)));
colSums(asvman);


################################################################################
#             2. Average the microbiota profiles for each host species 
#                             using a for loop --ALL HERBIVORES
################################################################################
#replace ASV table column names with their host species names
colnames(asvman)=meta$species_short[
  match(colnames(asvman),meta$Group)];

#average the microbiota profiles for each host species
myspec=as.character(levels(meta$species_short));
col_names=colnames(asvman); 

for(i in 1:length(myspec))
{
  A=grep(pattern=myspec[i], x=col_names);
  B=rowMeans(asvman[,A]);
  assign(paste(myspec[i]), B); 
};

#make master data frame of the averaged microbiota profiles
profiles=cbind(Cattle, Impala, Elephant, Zebra, Eland, Buffalo, 
               Giraffe,Warthog);

tprofiles=t(profiles); 
rowSums(tprofiles);


################################################################################
#             3. Average the microbiota profiles for each species 
#                             using a for loop -- BOVIDS ONLY
################################################################################
#simply remove non-bovid species from master data frame from above
tbovids=tprofiles[c(1,2,5,6),];


################################################################################
#             4. Calculate distance matrices for mantel test 
#                 ALL HERBIVORES
################################################################################
####BRAY-CURTIS distance
LMmantelbrayall=vegdist(tprofiles, method="bray"); 

####JACCARD distance
j=(tprofiles>0)*1
LMmanteljacall=vegdist(j, method="jaccard");

###WEIGHTED UNIFRAC distance
#first load phylogenetic tree of ASV sequences needed for Unifrac distances
load("data/Laikipia/04_LM_ASV_phylotree.Rdata");
uni<-phyloseq(otu_table(tprofiles, taxa_are_rows=FALSE),
              phy_tree(fitGTR$tree));
LMmantelunifracall=UniFrac(uni, 
                           weighted=TRUE, 
                           normalized=TRUE);

###UNWEIGHTED UNIFRAC distance
LMmantelunwunifracall=UniFrac(uni, 
                              weighted=FALSE, 
                              normalized=TRUE);


################################################################################
#             5. Calculate distance matrices for mantel test 
#                 BOVIDS ONLY
################################################################################

####BRAY-CURTIS distance
LMmantelbraybov=vegdist(tbovids, method="bray"); 

####JACCARD distance
j=(tbovids>0)*1
LMmanteljacbov=vegdist(j, method="jaccard");

###WEIGHTED UNIFRAC distance
uni<-phyloseq(otu_table(tbovids, taxa_are_rows=FALSE),
              phy_tree(fitGTR$tree));
LMmantelunifracbov=UniFrac(uni, 
                           weighted=TRUE, 
                           normalized=TRUE);

###UNWEIGHTED UNIFRAC distance
LMmantelunwunifracbov=UniFrac(uni, 
                              weighted=FALSE, 
                              normalized=TRUE);


################################################################################
#             6. Save distance matrices for mantel test 
################################################################################

save(LMmantelbrayall, LMmantelbraybov, LMmanteljacall, LMmanteljacbov, LMmantelunifracall, 
     LMmantelunifracbov, LMmantelunwunifracall, LMmantelunwunifracbov, 
     file = "data/Laikipia/05_LM_distances_mantel.Rdata");
