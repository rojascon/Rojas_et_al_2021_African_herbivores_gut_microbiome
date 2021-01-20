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

##CODE FOR: for 
#A) identifying ASVs present > 90% of samples 
#B) heatmap of the relative abundances of ^^these ASVs
#C) identifying ASVs biased towards a particular host species

source(file="scripts/00_background.R"); #load necessary packages and specifications

################################################################################
#             1. Load filtered ASV abundance table, ASV taxonomy table and 
#                         filtered metadata table                 
################################################################################
tax=read.csv("data/00_ASV_taxonomy_silva_v132.csv", header=T, row.names=1);
load("data/01_ASV_table_filtered.Rdata");
load("data/01_sample_metadata_filtered.Rdata");

#attach taxonomy to the ASV table 
tax=tax%>% select("asv_genus");
asv_tax=as.data.frame(t(asvdf));
asv_tax=merge(asv_tax,tax,by="row.names"); 
asv_tax$Row.names=NULL;


################################################################################
#             2. Identify ASVs present > 90% of samples                 
################################################################################
#rename last column
colnames(asv_tax)[ncol(asv_tax)]="taxa";

#convert ASV abundance table to ASV presence/absence table
cp=asv_tax; 
rownames(cp)=cp$taxa;
cp$taxa=NULL;
cp=(cp>0)*1; 
cp=as.data.frame(cp);

#identify ASVs present indiscriminantly across 90% of samples
cutoff=round(0.90*ncol(cp));
cp1=cp[rowSums(cp)>cutoff,];
print(rownames(cp1));


################################################################################
#             3. Make a heatmap of these "core" ASVs present > 90% of samples                 
################################################################################
#calculate ASV relative abundances (instead of presence/absence)
heat=asv_tax;
rownames(heat)=heat$taxa; heat$taxa=NULL;
heat=as.matrix(heat);
heat=prop.table(heat,2);
print(colSums(heat));

#subset dataframe to only include the ASVs you identified above [the > 90% ones]      
top=as.character(rownames(cp1));
heat2=heat[rownames(heat) %in% top,];

#reorganize columns of ASV table for heatmap
#so that all buffalo samples are together, then all of the cattle samples, etc
ideal.order=meta[
  with(meta, order(species_short, sample_month)),
  ];

heat2=heat2[,as.character(ideal.order$Group)];

#create a metadata file specific for heatmap [this is required for package 'pheatmap']
ideal.order=ideal.order[,c(1,5)];
rownames(ideal.order)=ideal.order$Group;
ideal.order$Group=NULL;

#select color-palette
species_colors <- list(species_short = c("#8dd3c7","#ffffb3","#bebada","#fb8072",
                                         "#80b1d3","#fdb462","#b3de69","#fccde5",
                                         "#d9d9d9","#bc80bd","#ccebc5"));

heat_color<-colorRampPalette(c("cornsilk","maroon4"), space = "rgb")(nrow(heat2));

names(species_colors$species_short) <- unique(ideal.order$species_short);

#plot heatmap
pheatmap(
  mat               = heat2,
  color             = heat_color,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  annotation_col    = ideal.order,
  annotation_colors = species_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = ""
);

#save manually by first manually adjusting R Plots window to desired size
#then going to Plots > Export > Save as PDF
#save file in directory "figures", with filename= 04_heatmap_core_ASVs.pdf


################################################################################
#             4. Identify ASVs biased towards a particular species                 
################################################################################
#criteria: ASVs must be present >75% of samples for that host species and 
#absent in 97% of other samples (aka "present" in <3% of other samples)
#these ASVs are "biased" towards that host speices

#create for loop to run this code for the 11 herbivore species
myspec=as.character(levels(meta$species_short));

for(i in 1:length(myspec))
{
  print(paste("Identifying % ASVs biased towards:", myspec[i]));
  #pick your samples
  mysamples=as.character(meta$Group
                         [meta$species_short==myspec[i]]);      
  A=cp[,which(colnames(cp) %in% mysamples)];
  lim=round(0.75*ncol(A))
  #find ASVs present in >75% of those samples
  ASVs_A=A[rowSums(A)>=lim,];                                   
  B=cp[,-which(colnames(cp) %in% mysamples)];
  otlim=round(0.03*ncol(B));
  #find ASVs present in <3% of the other samples
  ASVs_B=B[rowSums(B)<=otlim,];                                
  candidates=merge(ASVs_A, ASVs_B, by=0);
  #calculate the % ASVs biased towards host species
  val=(nrow(candidates)/nrow(ASVs_A))*100; print(val);
};

