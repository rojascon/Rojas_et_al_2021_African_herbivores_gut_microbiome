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

##CODE FOR: 
#A) subsampling study samples to 17,000 sequences each using mothur software
#B) calculating 4 metrics of gut microbiota alpha-diversity using phyloseq 
#       ASV richness, Chao1 richness, Shannon diversity, Simpson's diversity
#C) calculating Faith's PD using picante

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Subsample samples to 17000 sequence using mothur   
#               this number was chosen because it was the 2nd lowest number of 
#                             sequences found in a sample 
################################################################################

#load ASV table and ASV taxonomy
load("data/01_ASV_table_filtered.Rdata");
tax=read.csv("data/00_ASV_taxonomy_silva_v132.csv", header=T, row.names=1);

#make a df of ASV sequences and ASV labels
names=data.frame("seqs" = colnames(asvdf), 
                 "names" = sprintf("ASV%s",
                                   seq(1:ncol(asvdf))))
                                   
#format ASV table for mothur
#step1: replace colnames of ASV table with ASV labels from above
#step2:add the label, Group, and numOTUs columns to ASV table
colnames(asvdf)=names$names[match(names(asvdf),names$seqs)];
tempdf=data.frame("label" = rep(0.03,nrow(asvdf)), 
                  "Group" = row.names(asvdf), 
                  "numOtus" = rep(ncol(asvdf),nrow(asvdf)))
mothur=cbind(tempdf, asvdf)
write.table(mothur, file="data/07_input.mothur.txt", row.names=F, sep="\t")

#run mothur
#step1: for some reason, must open the above .txt file in Excel, and 
#           save as Tab delimited file (.txt)
#step2: download and open mothur software (Schloss et.al 2009)
#step3: run the command: 
#           sub.sample(shared=07_input.mothur.txt, size=17000, persample=true)
#step4: rename output file as "07_output.mothur.txt") and place in the data directory


################################################################################
#             2. Clean file output from mothur
################################################################################

#read in output from mothur
mothur2=read.table("data/07_output.mothur.txt", sep="\t",header=T);
rownames(mothur2)=mothur2$Group
mothur2=mothur2[,-c(1:3)]

#replace ASV labels with ASV sequence names
colnames(mothur2)=names$seqs[match
                             (names(mothur2),names$names)];
mothur2=as.matrix(mothur2);


################################################################################
#             3. Calculate gut microbiota alpha-diversity using phyloseq
################################################################################
#make a phyloseq object and calculate 4 metrics of alpha-diversity
ps<- phyloseq(otu_table(mothur2, taxa_are_rows=FALSE));
alpha=estimate_richness(ps,split = TRUE, 
                        measures = c("Observed","Chao1","Shannon","Simpson"));

#calculate Good's coverage and append to alpha df
gcov=goods(mothur2);
alpha=merge(alpha, gcov, by="row.names", all=TRUE); 
colnames(alpha)[1]="Group";


################################################################################
#             4. Calculate Faith's Phylogenetic diversity using picante
################################################################################

#load phylogenetic tree of ASVs sequences
#see script "05.1_betadiversity_distances.R" for how to generate
load("data/05_ASV_phylotree.Rdata"); 

#calculate Faith's phylogenetic diversity and append to alpha df
faith=pd(mothur2, phy_tree(fitGTR$tree), include.root=F);
faith$Group=row.names(faith);
alpha=inner_join(alpha, faith, by="Group");

#save data frame of alpha diversity metrics 
save(alpha, file="data/07_alpha_values.Rdata");
