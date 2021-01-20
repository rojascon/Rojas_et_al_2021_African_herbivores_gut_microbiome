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

##CODE FOR: filtering ASV abundance table and metadata to only include
##samples of interest

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load ASV abundance table, ASV taxonomy table and 
#                         metadata table                 
################################################################################
#load ASV table output by DADA2 R tutorial 
#https://benjjneb.github.io/dada2/tutorial.html
load("data/Laikipia/00_LM_seqtab_dada2.Rdata");

#load metadata and ASV taxonomy
meta=read.csv("data/Laikipia/00_LM_sample_metadata.csv", stringsAsFactors = F);
tax=read.csv("data/Laikipia/00_LM_ASV_taxonomy_silva_v132.csv", header=T, row.names=1); 


################################################################################
#             2. Remove Baboon samples and samples that did not sequence well 
#                           from metadata table
################################################################################
meta=meta[meta$species_short!="Baboon",];

bad_samples=c("AN102", "AN11","AN120","AN127","AN128","AN159","AN166","AN191",
"AN194","AN198","AN27","AN44","AN45");

meta=meta[!meta$Group %in% bad_samples,];

#reorder factors of metadata
meta$sample_month=factor(meta$sample_month, levels=c("march","april","may","june","july",
                                                         "august","october","november"));
meta$species_short=factor(meta$species_short, levels=c("Buffalo","Cattle","Eland","Impala",
                                                           "Giraffe","Warthog","Zebra","Elephant"));
meta$diet_guild=factor(meta$diet_guild, 
                         levels=c("grazer","browser","mixed_feeder"));
meta$region=factor(meta$region, 
                     levels=c("Masai_Mara","Laikipia"));
meta=meta[order(meta$Group),];


################################################################################
#             3. Remove Baboon samples and samples that did not sequence well 
#                             from ASV table
################################################################################
#convert ASV table from matrix to data frame
asvdf=as.data.frame(seqtab.nochim);

#remove baboon samples and samples that did not sequence well
asvdf=asvdf[rownames(asvdf) %in% meta$Group,];

#remove ASVs classified as Eukarya, Chloroplast, Mitochondria, and Unknown
#these have already been removed from the taxonomy file
asvdf=asvdf[,colnames(asvdf) %in% rownames(tax)];

#remove ASVs that are not present in any samples
asvpa=(t(asvdf));
asvpa=(asvpa>0)*1; 
asvpa=as.data.frame(asvpa);
asvpa$sm=rowSums(asvpa);
asvpa=asvpa[asvpa$sm>0,];

asvdf=asvdf[,colnames(asvdf) %in% rownames(asvpa)];


################################################################################
#             4. save filtered ASV table and filtered meta_data
################################################################################
save(asvdf, file="data/Laikipia/01_LM_ASV_table_filtered.Rdata");
save(meta, file="data/Laikipia/01_LM_sample_metadata_filtered.Rdata");

