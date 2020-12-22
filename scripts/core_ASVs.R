#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for identifying ASVs present > 90% of samples 
#Code for identifying ASVs present only in one host species

source(file="scripts/background.R"); #load necessary packages and specifications

#load metadata with baboon samples removed, as well as samples that did not amplify
load("data/sample_metadata_beta.Rdata");

#load ASV table and ASV taxonomy output by DADA2 
#samples that did not amplify were removed; see script "get_ASV_table.R"
#chloroplast, mitochondria, and unknown were removed from taxonomy
load("data/ASV_table.Rdata");
tax=read.csv("data/ASV_taxonomy_silva_v132.csv", header=T); 

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

######################  IDENTIFY ASVs present across 90% of samples  ###############
#trim unecessary columns from ASV taxonomy abundance table
core=asv_tax[,which(names(asv_tax) 
                    %in% c(samples, "asv_genus"))];
colnames(core)[ncol(core)]="taxa";

#convert ASV abundance table to ASV presence/absence
cp=core; rownames(cp)=cp$taxa; cp$taxa=NULL;
cp=(cp>0)*1; 
cp=as.data.frame(cp);

#identify ASVs present indiscriminantly across 90% of samples
cutoff=round(0.90*ncol(cp));
cp1=cp[rowSums(cp)>cutoff,];
print(rownames(cp1));

###################### HEATMAP of the relative abundance of these ASVs  ###############

#calculate ASV relative abundances (instead of presence/absence)
heat=core;
rownames(heat)=heat$taxa; heat$taxa=NULL;
heat=as.matrix(heat);
heat=prop.table(heat,2);
print(colSums(heat));

#subset dataframe to only include the ASVs you identified above [the > 90% ones]      
top=as.character(rownames(cp1));
heat2=heat[rownames(heat) %in% top,];

#reorganize columns of ASV table for heatmap
#so that all buffalo samples are together, then all of the cattle samples, etc
ideal.order=metadf[
  with(metadf, order(species_short, sample_month)),
  ];

heat2=heat2[,as.character(ideal.order$Group)];

#prepare metadata for heatmap (each sample name and the herbivore species it comes from)
ideal.order=ideal.order[,c(1,5)];
rownames(ideal.order)=ideal.order$Group;
ideal.order$Group=NULL;

#select color-palette
species_colors <- list(species_short = c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3",
                                         "#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd",
                                         "#ccebc5"));

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
#set directory as "figures"

######################  IDENTIFY ASVs BIASED TOWARDS PARTICULAR HOST SPECIES ###############
#criteria: ASVs present >75% of samples for that host sp AND in <3% of the other samples
#these ASVs mostly associate with one particular host species

#find the ASVs present >75% of samples for a host species
mysamples=as.character(metadf$Group
                       [metadf$species_short=="Zebra"]); #select your host species
cpsp=cp[,which(colnames(cp) %in% mysamples)];
lim=round(0.75*ncol(cpsp));
myASVs=cpsp[rowSums(cpsp)>=lim,];

#find the ASVs present in <3% of remaining samples
other=as.character(metadf$Group
                   [metadf$species_short!="Zebra"]); #same host species here
cpot=cp[,colnames(cp) %in% other,];
otlim=round(0.03*ncol(cpot));
rareASVs=cpot[rowSums(cpot)<=otlim,];

#find the overlap in the two data frames
# aka find ASVs present >75% of samples for that host sp 
#AND in <3% of remaining samples
list_of_data = list(myASVs,rareASVs);
finalASVs= Reduce(intersect, 
                  lapply(list_of_data, row.names)); 

#calculate the percent of ASVs biased towards that host species
#basically dividing the number of ASVs that meet both criteria / by the number of ASVs that meet the first criterion
value1=length(finalASVs)/
  nrow(myASVs); 

print(value1*100);
