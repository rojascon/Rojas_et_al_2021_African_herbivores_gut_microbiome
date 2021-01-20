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
#working from Laikipia subdirectories within the main scripts, 
#data, and figures directories

##CODE FOR: 
#A) identifying ASVs present > 90% of samples
#B) generating heatmap of 32 most abundant ASVs found across samples
#C) running linear models to statistically confirm patterns from heatmap

source(file="scripts/00_background.R"); #load necessary packages and specifications

################################################################################
#             1. Load filtered ASV abundance table, ASV taxonomy table and 
#                         filtered metadata table                 
################################################################################
tax=read.csv("data/Laikipia/00_LM_ASV_taxonomy_silva_v132.csv", header=T, row.names=1);
load("data/Laikipia/01_LM_ASV_table_filtered.Rdata");
load("data/Laikipia/01_LM_sample_metadata_filtered.Rdata");

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
#             3. Make a heatmap of the 32 most abundant ASVs                 
################################################################################
#calculate ASV relative abundances (instead of presence/absence)
heat=asv_tax;
rownames(heat)=heat$taxa; heat$taxa=NULL;
heat=as.matrix(heat);
heat=prop.table(heat,2);
#print(colSums(heat));

#keep the top 32 most abundant ASVs (or ASVs with >0.23% relative abundance)
heat2=as.data.frame(heat);
heat2$AVG=rowMeans(heat2);
heat2=heat2[heat2$AVG>0.0023,];
heat2=heat2[order(-heat2$AVG),];
heat2$AVG=NULL; 
print(paste("the number of ASVs that will show on heatmap is:", nrow(heat2), sep=" "));

#reorganize columns of ASV table for heatmap
#so that masai mara buffalo, are followed by laikipia buffalo, then masai mara cattle,etc
ideal.order=meta[
  with(meta, order(species_short, region)),
  ];

heat2=heat2[,as.character(ideal.order$Group)];

#create a metadata file specific for heatmap [this is required for package 'pheatmap']
ideal.order=ideal.order[,c(1,15)];
rownames(ideal.order)=ideal.order$Group;
ideal.order$Group=NULL;

#select color-palettes
region_colors <- list(region = c("#a1d99b","#beaed4"));

heat_color<-colorRampPalette(c("cornsilk","maroon4"), 
                             space = "rgb")(nrow(heat2));

names(region_colors$region) <- unique(ideal.order$region);


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
  annotation_colors = region_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = ""
);

#save manually by first manually adjusting R Plots window to desired size
#then going to Plots > Export > Save as PDF
#save file in directory "figures/Laikipia", with filename= 02_laikipia_mara_ASVheatmap.pdf

################################################################################
#             3. Run Linear Mixed Models to test whether ASV abundances
#                   shown in the heatmap differ between hosts from the two 
#                           geographic regions
################################################################################
# formula will be: 
#ASV[i] ~ host geographic region + (1|host species)
#restricting analyses to WITHIN a host species

#prepare dataframe for analyses
asvtest=as.data.frame(t(heat2));
rownames(meta)=meta$Group; meta$Group=NULL;
asvtest=merge(asvtest, meta, by="row.names");
rownames(asvtest)=asvtest$Row.names; 
asvtest$Row.names=NULL;

### run models using a for-loop [32 ASVs to test, so 32 models]
#but first set up this empty data frame to store model output
mypvalues <- data.frame(index = seq_len(nrow(x = heat2)+1),
                        ASV=NA,
                        estimate=NA,
                        stderror=NA,
                        pvalue = NA);

for(i in 1:nrow(heat2)) 
{
  asvmodel=lmer(asvtest[,i]~region + (1|species_short), data=asvtest);
  estimate=coef(summary(asvmodel))[,1];
  sterror=coef(summary(asvmodel))[,2];
  pval=coef(summary(asvmodel))[,5];
  mypvalues[i+1, 2] <-colnames(asvtest)[i];
  mypvalues[i+1, 3] <-estimate[2];
  mypvalues[i+1, 4] <-sterror[2];
  mypvalues[i+1, 5] <-pval[2];
};

#add p-values to your table and adjust for multiple comparisons
mypvalues=mypvalues[-1,-1];
mypvalues$pvalue_adjusted=p.adjust(mypvalues$pvalue, method="BH");

#export table
write.csv(mypvalues, file="figures/Laikipia/02_LM_lmer_ASVs.csv", row.names=F);