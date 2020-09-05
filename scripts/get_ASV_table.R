#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for finalizing ASV table and removing samples that did not amplify

source(file="scripts/background.R"); #load necessary packages and specifications

#load ASV table output by DADA2 R tutorial & file that contains the # reads/sample
#https://benjjneb.github.io/dada2/tutorial.html
load("data/seqtab_dada2.Rdata");
reads=read.csv("data/numseqs_dada2.csv", header=T);

#clean up sample names in ASV table
asv.tbl=as.data.frame(seqtab.nochim);
rownames(asv.tbl)=gsub("_.*", "", rownames(asv.tbl));

#remove samples that did not amplify well (<400 reads)
mysamples=reads$Group[reads$nonchim>400];
asv.tbl=asv.tbl[which(row.names(asv.tbl) %in% mysamples),];
save(asv.tbl, file="data/ASV_table.Rdata");