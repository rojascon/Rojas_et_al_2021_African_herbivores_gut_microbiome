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

##CODE FOR: identifying the bacterial families that are significantly 
#associated with particular:
#A) herbivore families
#B) herbivore dietary guilds 
# as determined by indicator species analysis

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load filtered ASV abundance table, ASV taxonomy table and 
#                         filtered metadata table                 
################################################################################
tax=read.csv("data/00_ASV_taxonomy_silva_v132.csv", header=T, row.names=1);
load("data/01_ASV_table_filtered.Rdata");
load("data/01_sample_metadata_filtered.Rdata");

#attach taxonomy to the ASV table 
tax=tax%>% select("Family");
asv_tax=as.data.frame(t(asvdf));
asv_tax=merge(asv_tax,tax,by="row.names"); 
asv_tax$Row.names=NULL;


################################################################################
#             2. Calculate bacterial family relative abundances before
#               running indicator species analysis
################################################################################

#calculate bacterial family relative abundances 
fam=aggregate(.~Family, asv_tax, sum);  
fam[,-1] <- lapply(fam[,-1], function(x) (x/sum(x))*100);
print(colSums(fam[-1]));

#transpose table of bacterial family abundances
rownames(fam)=fam$Family;
fam$Family=NULL;
famt=as.data.frame(t(fam));


################################################################################
#             3. Determine bacterial families associated with particular
#                 herbivore families using indicator species analysis
################################################################################

#get list of bacterial taxa that may be differentially abundant among host families
#pay attention only to results from: Group Giraffidae, Group Suidae, 
#Group Equidae, Group Elephantidae (NO INTERACTION TERMS)

indbac=multipatt(famt, meta$Family, control=how(nperm=999));
#print(summary(indbac));

#save findings to a text file
sink("data/03_indicator_taxa_host_families.txt");
summary(indbac);
sink(); 


################################################################################
#             4. Determine bacterial families associated with particular
#                 herbivore dietary guilds using indicator species analysis
################################################################################

#get list of bacterial taxa that may be differentially abundant among host dietary guilds
#pay attention only to results from: Group grazers, Group browsers, Group mixed_feeders
#NO INTERACTION TERMS

indbac2=multipatt(famt, meta$diet_guild, control=how(nperm=999));
#print(summary(indbac2));

#save findings to a text file
sink("data/03_indicator_taxa_host_dietguilds.txt");
summary(indbac2);
sink(); 
