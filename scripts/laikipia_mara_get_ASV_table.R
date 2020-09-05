#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#MASAI MARA && LAIKIPIA DATA COMBINED
#Code for finalizing ASV table/metadata and removing samples that did not amplify/are baboon samples

source(file="scripts/background.R"); #load necessary packages and specifications

#load ASV table output by DADA2 R tutorial & file that contains the # reads/sample
#https://benjjneb.github.io/dada2/tutorial.html
load("data/laikipia_mara_seqtab_dada2.Rdata");
reads=read.csv("data/laikipia_mara_numseqs_dada2.csv", header=T);
colnames(reads)[1]="Group";

#load combined metadata file & remove baboon samples
metadf=read.csv("data/laikipia_mara_metadata.csv");
metadf=metadf[metadf$species_short!="Baboon",];

#remove samples that did not amplify well (<400 reads) from metadata & ASV table
mysamples=reads$Group[reads$nonchim>405];
metadf=metadf[metadf$Group %in% mysamples,];
asv.tbl=seqtab.nochim[which(row.names(seqtab.nochim) %in% metadf$Group),];

#reorder factors of metadata
metadf$sample_month=factor(metadf$sample_month, levels=c("march","april","may","june","july",
                                                         "august","october","november"));
metadf$species_short=factor(metadf$species_short, levels=c("Buffalo","Cattle","Eland","Impala",
                                                       "Giraffe","Warthog","Zebra","Elephant"));
metadf$diet_guild=factor(metadf$diet_guild, 
                         levels=c("grazer","browser","mixed_feeder"));
metadf$region=factor(metadf$region, 
                     levels=c("Masai_Mara","Laikipia"));
metadf=metadf[order(metadf$Group),];

#save formatted ASV table and metadata
save(asv.tbl, file="data/laikipia_mara_filtered_ASV_table.Rdata");
save(metadf, file="data/laikipia_mara_filtered_metadata.Rdata");

