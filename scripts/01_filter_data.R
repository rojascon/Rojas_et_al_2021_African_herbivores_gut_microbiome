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

##CODE FOR: filtering ASV abundance table and metadata to only include
##samples of interest

source(file="scripts/00_background.R"); #load necessary packages and specifications

################################################################################
#             1. Load ASV abundance table, ASV taxonomy table and 
#                         metadata table                 
################################################################################
#load ASV table output by DADA2 R tutorial 
#https://benjjneb.github.io/dada2/tutorial.html
load("data/00_seqtab_dada2.Rdata");

#load metadata and ASV taxonomy
meta=read.csv("data/00_sample_metadata.csv", stringsAsFactors = F);
tax=read.csv("data/00_ASV_taxonomy_silva_v132.csv", header=T, row.names=1); 

################################################################################
#             2. Remove Baboon samples and samples that did not sequence well 
#                 from metadata table; format metadata factors
################################################################################
meta=meta[meta$species_short!="Baboon",];

bad_samples=c("AN159","AN194","AN60", "AN83","AN118", "AN154", "AN36",  "AN141", 
              "AN102", "AN127", "AN45","AN11","AN198", "AN128","AN27","AN120");

meta=meta[!meta$Group %in% bad_samples,];

#format metada file factors
meta$sample_month=factor(meta$sample_month, 
                         levels=c("march","april","may","june"));

meta$Order=factor(meta$Order, 
                  levels=c("Artiodactyla","Perissodactyla","Proboscidea"));

meta$Family=factor(meta$Family, 
                   levels=c("Bovidae","Giraffidae","Suidae","Equidae","Elephantidae"));

meta$species_short=factor(meta$species_short, 
                          levels=c("Buffalo","Cattle","Topi","Eland","Impala","Gazelle",
                                   "Dikdik","Giraffe","Zebra","Warthog","Elephant"));

meta$diet_guild=factor(meta$diet_guild, 
                       levels=c("grazer","browser","mixed_feeder"));

meta=meta[order(meta$Group),];


################################################################################
#             3. Remove Baboon samples and samples that did not sequence well 
#                 from ASV table;
################################################################################
#clean up sample names in ASV table
asvdf=as.data.frame(seqtab.nochim);
rownames(asvdf)=gsub("_.*", "", rownames(asvdf));

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
save(asvdf, file="data/01_ASV_table_filtered.Rdata");
save(meta, file="data/01_sample_metadata_filtered.Rdata");
