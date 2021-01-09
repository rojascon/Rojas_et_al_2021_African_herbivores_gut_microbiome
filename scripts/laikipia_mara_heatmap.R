#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#MASAI MARA && LAIKIPIA DATA COMBINED

#Code for generating heatmap of top 30 most abundant ASVs across samples
##distinguishing between samples from Masai Mara and. Laikipia
##also code for running linear models to statistically confirm patterns from heatmap

source(file="scripts/background.R"); #load necessary packages and specifications

#load ASV table output by DADA2 && metadata
#samples that did not amplify & baboon samples were removed
#see script "laikipia_mara_get_ASV_table.R"
load("data/laikipia_mara_filtered_ASV_table.Rdata");
load("data/laikipia_mara_filtered_metadata.Rdata");
nrow(metadf)==nrow(asv.tbl);

#load ASV taxonomy output by DADA2
tax=read.csv("data/laikipia_mara_ASV_taxonomy_silva_v132.csv", header=T); 

#make a combined ASV & genus label (e.g. ASV1|Treponema)
tax$label=paste("ASV",sep="",rownames(tax));
tax$asv_genus=paste(tax$label,sep="|",tax$Genus);
rownames(tax)=tax$ASV; 
tax$ASV=NULL;

#attach taxonomy to the ASV table 
asvt=as.data.frame(t(asv.tbl));
asv_tax=merge(asvt,tax,by="row.names"); 
colnames(asv_tax)[1]="ASV";

#make a vector of the samples you want for analysis
samples=as.character(metadf$Group);

#######################  IDENTIFY ASVs present across 90% of samples ##################
#trim unecessary columns from ASV abundance table
core=asv_tax[,which(names(asv_tax) 
                    %in% c(samples, "asv_genus"))];
colnames(core)[ncol(core)]="taxa";

#convert ASV abundance table to ASV presence/absence
cp=core; rownames(cp)=cp$taxa; cp$taxa=NULL;
cp$seqs=rowSums(cp); cp=cp[cp$seqs>2,]; cp$seqs=NULL;
cp=(cp>0)*1; 
cp=as.data.frame(cp);

#identify ASVs present indiscriminantly across 90% of samples
cutoff=round(0.90*ncol(cp));
cp1=cp[rowSums(cp)>cutoff,];

print(rownames(cp1));

#######################  Heatmap of 32 most abundant ASVs  #######################
#trim unecessary columns from ASV taxonomy abundance table
asvh=asv_tax[,which(names(asv_tax) 
                    %in% c(samples, "asv_genus"))];
colnames(asvh)[ncol(asvh)]="taxa";
rownames(asvh)=asvh$taxa; asvh$taxa=NULL;
asvh=as.matrix(asvh);

#calculate ASV relative abundances 
asvh=prop.table(asvh,2);
print(colSums(asvh));

#keep the top 32 most abundant ASVs (or ASVs with >0.23% relative abundance)
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
region_colors <- list(region = c("#a1d99b","#beaed4"));

heat_color<-colorRampPalette(c("cornsilk","maroon4"), 
                             space = "rgb")(nrow(asv_heat));

names(region_colors$region) <- unique(ideal.order$region);

#plot heatmap
pheatmap(
  mat               = asv_heat,
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

##save heatmap
#dev.off();

#save manually by going to Plots > Export > Save as PDF

#######################  LINEAR MODELS TO TEST 32 MOST ABUNDANT ASVs  ##############################
################### example: ASV 10 ~ host geographic region ##############

##only testing the top 32 abundant ASVs from above
#same ASVs from the heatmap

#prepare dataframe for analyses
asvtest=as.data.frame(t(asv_heat));
rownames(metadf)=metadf$Group; metadf$Group=NULL;
asvtest=merge(asvtest, metadf, by=0);
rownames(asvtest)=asvtest$Row.names; 
asvtest$Row.names=NULL;

### run linear mixed effects model on all ASVs at the same time using a for-loop
## MODEL: y (ASV 12315) ~ host region + (1|host species), data=asvtest

#set up this empty data frame to store model pvalues
mypvalues <- data.frame(index = seq_len(nrow(x = asv_heat)+1),
                        ASV=NA,
                        estimate=NA,
                        stderror=NA,
                        pvalue = NA);

for(i in 1:nrow(asv_heat)) 
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
write.csv(mypvalues, file="figures/lmer_topASVs_hostregion.csv", row.names=F);
