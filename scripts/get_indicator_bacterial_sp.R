#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for identifying the list of bacterial taxa that are enriched in particular groups
#using indicator species analysis

source(file="scripts/background.R"); #load necessary packages and specifications

#load metadata with baboon samples removed, as well as samples that did not amplify
load("data/sample_metadata_beta.Rdata");

#load ASV taxonomy output by DADA2 
#chloroplast, mitochondria, and unknown were removed 
tax=read.csv("data/ASV_taxonomy_silva_v132.csv", header=T); 
rownames(tax)=tax$ASV; tax$ASV=NULL;

#load ASV table output by DADA2 
#samples that did not amplify were removed; see script "get_ASV_table.R"
load("data/ASV_table.Rdata");

#remove baboon samples from ASV table
asvdf=asv.tbl[rownames(asv.tbl) %in% metadf$Group,]; 
nrow(metadf)==nrow(asvdf);

#attach taxonomy to the ASV table 
asv_tax=as.data.frame(t(asvdf));
asv_tax=merge(asv_tax,tax,by="row.names"); 
colnames(asv_tax)[1]="ASV";


#######################  PREPARE ASV TABLE FOR INDICATOR SPECIES ANALYSIS ##############################
################# TO IDENTIFY BACTERIAL FAMILIES THAT MAY BE ENRICHED IN PARTICULAR HOST GROUPS #################

#select bacterial taxonomic rank 
family=asv_tax[,which(names(asv_tax) 
                      %in% c(as.character(metadf$Group), "Family"))];
colnames(family)[ncol(family)]="taxa";

#calculate bacterial family relative abundances 
fam=aggregate(.~taxa, family, sum);  
fam[,-1] <- lapply(fam[,-1], function(x) (x/sum(x))*100);
print(colSums(fam[-1]));

#transpose table of bacterial family abundances
rownames(fam)=fam$taxa;
fam$taxa=NULL;
famt=as.data.frame(t(fam));
nrow(famt)==nrow(metadf)

####################### RUN INDICATOR SPECIES ANALYSIS ##############################

#get list of bacterial taxa that may be differentially abundant among host families
#pay attention only to results from: Group Giraffidae, Group Suidae, 
#Group Equidae, Group Elephantidae
indbac=multipatt(famt, metadf$Family, control=how(nperm=999))
summary(indbac)

#save findings to a text file
sink("data/indicator_taxa_host_families.txt");
print(summary(indbac));
sink(); 

#get list of bacterial taxa that may be differentially abundant among host dietary guilds
#pay attention only to results from: Group grazers, Group browsers, Group mixed_feeders
indbac=multipatt(famt, metadf$diet_guild, control=how(nperm=999))
summary(indbac)

#save findings to a text file
sink("data/indicator_taxa_host_dietguilds.txt");
print(summary(indbac));
sink(); 




