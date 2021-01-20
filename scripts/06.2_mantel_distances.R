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

##CODE FOR: generating microbiota distance matrices Bray-Curtis (bray), 
#Jaccard (jac), Weighted Unifrac (wuni), Unweighted Unifrac (unwuni) for Mantel tests

#These distances are based on microbiota profiles that were AVERAGED across host species 
#Their dimensions match the dimensions of the matrix of host divergence times

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load filtered ASV table & filtered metadata table                 
################################################################################
load("data/01_sample_metadata_filtered.Rdata");
load("data/01_ASV_table_filtered.Rdata");

#turn ASV counts to ASV proportions
asvman<- apply(asvdf, 1, function(i) (i/sum(i)));
colSums(asvman);


################################################################################
#             2. Average the microbiota profiles for each species 
#                       using a for loop -- ALL HERBIVORES
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
profiles=cbind(Cattle, Dikdik, Impala, Topi, Elephant, Zebra, Eland, Buffalo, 
               Giraffe,Gazelle, Warthog);

tprofiles=t(profiles); 
rowSums(tprofiles);

################################################################################
#             3. Average the microbiota profiles for each species 
#                       using a for loop -- BOVIDS ONLY
################################################################################
#simply remove non-bovid species from master data frame from above
tbovids=tprofiles[c(1:4,7:8,10),];


################################################################################
#             4. Calculate distance matrices for mantel test --ALL HERBIVORES
################################################################################
####BRAY-CURTIS distances
mantelbrayall=vegdist(tprofiles, method="bray"); 

####JACCARD distances
j=(tprofiles>0)*1
manteljacall=vegdist(j, method="jaccard");

###WEIGHTED UNIFRAC distances
#first load phylogenetic tree of ASV sequences needed for Unifrac distances
load("data/05_ASV_phylotree.Rdata");
uni<-phyloseq(otu_table(tprofiles, taxa_are_rows=FALSE),
              phy_tree(fitGTR$tree));
mantelunifracall=UniFrac(uni, 
                         weighted=TRUE, 
                         normalized=TRUE);

###UNWEIGHTED UNIFRAC distances
mantelunwunifracall=UniFrac(uni, 
                            weighted=FALSE, 
                            normalized=TRUE);


################################################################################
#             5. Calculate distance matrices for mantel test --BOVIDS ONLY
################################################################################
####BRAY-CURTIS distances
mantelbraybov=vegdist(tbovids, method="bray"); 

####JACCARD distances
j=(tbovids>0)*1;
manteljacbov=vegdist(j, method="jaccard");

###WEIGHTED UNIFRAC distances
uni<-phyloseq(otu_table(tbovids, taxa_are_rows=FALSE),
              phy_tree(fitGTR$tree));
mantelunifracbov=UniFrac(uni, 
                         weighted=TRUE, 
                         normalized=TRUE);

###UNWEIGHTED UNIFRAC distances
mantelunwunifracbov=UniFrac(uni, 
                            weighted=FALSE, 
                            normalized=TRUE);


################################################################################
#             6. Save distance matrices for mantel test 
################################################################################
save(mantelbrayall, mantelbraybov, manteljacall, manteljacbov, mantelunifracall, 
     mantelunifracbov, mantelunwunifracall, mantelunwunifracbov, 
     file = "data/06_distances_mantel.Rdata");
