#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#MASAI MARA && LAIKIPIA DATA COMBINED
#Code for heatmap of gut microbiota compositions--top 30 most abundant ASVs

source(file="scripts/background.R"); #load necessary packages and specifications

#load ASV table output by DADA2 && metadata
#samples that did not amplify & baboon samples were removed
#see script "laikipia_mara_get_ASV_table.R"
load("data/laikipia_mara_filtered_ASV_table.Rdata");
load("data/laikipia_mara_filtered_metadata.Rdata");
nrow(metadf)==nrow(asv.tbl);

#load ASV taxonomy output by DADA2
tax=read.csv("data/laikipia_mara_ASV_taxonomy_silva_v132.csv", header=T); 

#add an ASV number column to taxonomy ("ASV1", "ASV2", etc)
#to this ASV label, attach the bacterial genus ("ASV1|")
tax_heat=tax; 
tax_heat$label=paste("ASV",sep="",rownames(tax_heat));
tax_heat$label_genus=paste(tax_heat$label,sep="|",tax_heat$Genus);
rownames(tax_heat)=tax_heat$ASV; 
tax_heat$ASV=NULL;

#attach taxonomy to the ASV table 
asvt=as.data.frame(t(asv.tbl));
asv_tax=merge(asvt,tax_heat,by="row.names"); 
colnames(asv_tax)[1]="ASV";

#make a vector of the samples you want for analysis
samples=as.character(metadf$Group);

#######################  ASV counts to ASV proportions  ##############################
#trim unecessary columns from ASV taxonomy abundance table
asvh=asv_tax[,which(names(asv_tax) 
                      %in% c(samples, "label_genus"))];
rownames(asvh)=asvh$label_genus;
asvh$label_genus=NULL;
asvh=as.matrix(asvh);

#calculate ASV relative abundances 
asvh=prop.table(asvh,2);
print(colSums(asvh));

#keep phyla >0.23% relative abundance across samples
##OR the top 30 most abundant ASVs
asv_heat=as.data.frame(asvh);
asv_heat$AVG=rowMeans(asv_heat);
asv_heat=asv_heat[asv_heat$AVG>0.0023,];
asv_heat=asv_heat[order(-asv_heat$AVG),];
asv_heat$AVG=NULL; 
paste("the number of ASVs that will show on heatmap is:", nrow(asv_heat), sep=" ");

#reorganize columns of ASV table for heatmap
#masai mara buffalo, laikipia buffalo, masai mara cattle, laikipia cattle, etc
ideal.order=metadf[
  with(metadf, order(species_short, region)),
  ];
asv_heat=asv_heat[,as.character(ideal.order$Group)];

#remove unecessary columns from metadata
ideal.order=ideal.order[,c(1,15)];
rownames(ideal.order)=ideal.order$Group;
ideal.order$Group=NULL;

#select color-palette
mat_colors <- list(region = c("#a1d99b","#beaed4"));
myCol<-colorRampPalette(c("cornsilk","maroon4"), space = "rgb")(nrow(asv_heat));
names(mat_colors$region) <- unique(ideal.order$region);

#prepare plot for saving
pdf(file="figures/laikipia_mara_ASVheatmap.pdf", width=16, height=8);

#plot heatmap
pheatmap(
  mat               = asv_heat,
  color             = myCol,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  annotation_col    = ideal.order,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = ""
);

##save heatmap
dev.off();




