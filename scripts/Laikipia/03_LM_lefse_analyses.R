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
#A) generating ASV tables for LEfSe analyses
#               LEfSe will identify the ASVs that are differentially enriched 
#                   in Masai Mara giraffes vs Laikipia giraffes, for example
#B) plotting the output from LEfSe in the form of diverging plots

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
rownames(asv_tax)=asv_tax$asv_genus
asv_tax$Row.names=NULL;
asv_tax$asv_genus=NULL;


################################################################################
#             2. Calculate ASV relative abundances                 
################################################################################

asvlef=as.matrix(asv_tax);
asvlef=prop.table(asvlef,2);
#print(colSums(asvlef));


################################################################################
#             3. Make an ASV table for each host species, so that only
#                 samples from that host species are included
################################################################################

#samples will be labelled by their host's geographic region
#only ASVs >0.01% average relative abundance for that species are included in data frame

myspecies=levels(meta$species_short);

for(i in 1:length(myspecies))
{
  A=meta$Group[meta$species_short==myspecies[i]];
  B=as.data.frame(asvlef[, colnames(asvlef) %in% A]);
  colnames(B)=meta$region[
    match(colnames(B),meta$Group)];
  B$avg=rowMeans(B);
  B=B[B$avg>0.0001,];
  B$avg=NULL;
  assign(paste(myspecies[i]), B); 
};

#save output files from for-loop as .txt for LEfSe
mydfs=list(Buffalo, Cattle, Eland, Impala, Giraffe, Warthog, Zebra, Elephant);

for (i in 1:length(myspecies)){
  filename <- paste("data/Laikipia/03_LM_",myspecies[i],"_lefse.txt", sep="");
  write.table(mydfs[[i]], filename, col.names=NA,row.names=TRUE,sep="\t")
};

################################################################################
#             3. Run LEFSE on the Galaxy platform
#               https://huttenhower.sph.harvard.edu/galaxy/
################################################################################

#Step1: Upload data frames from above by going to Get Data > 
# upload file from your Computer

#Step2: Go to LEfSe > Format Data for LEfSe 
#Select whether the vectors (features and meta-data information) are listed in 
#rows or columns: select ROWS
#Select which row to use as class: use #1 region

#Step3: set the Alpha value for the factorial Kruskal-Wallis test among classes = 0.05
#Step4: compile all tables output by LEFSE into one file = "03_LM_lefse_output.csv"
#Place that file in the data/Laikipia folder


################################################################################
#             3. Prepare LEFSE data for plotting inggplot2
################################################################################

#read in LefSe output
lefse=read.csv("data/Laikipia/03_LM_lefse_output.csv", header=T);
lefse[,7:9]=str_split_fixed(lefse$ASV, "_",3); lefse$V9=NULL;

colnames(lefse)=c("ASV","mean","region","LDA","pvalue","species","label","taxa");

lefse$species=factor(lefse$species, levels=c("Buffalo", "Cattle", "Eland", 
                                             "Impala", "Giraffe", "Warthog", 
                                             "Zebra", "Elephant"));

#adjust pvalues for multiple comparisons
lefse$padjusted=p.adjust(lefse$pvalue, method="BH");
lefse=lefse[lefse$label!="NA",];

#only retain ASVs with LDA effect sizes > 3.2 for plotting purposes
#you can set this cutoff as any value between 2 and 4.30
lefse2=lefse[lefse$padjusted<0.05,];
lefse2=lefse2[lefse2$LDA>3.2,];

#convert all Laikipia effect sizes to -values for diverging plots
lefse2$LDA[lefse2$region=="Laikipia"]=
  (lefse2$LDA[lefse2$region=="Laikipia"])*-1

################################################################################
#             3. Plot output from LEFSE as diverging plots
#               to show which ASVs are differentially enriched in 
#                   Masai Mara relative to Laikipia herbivores
################################################################################

#set up color palette and LDA tickmark labels for plot
region_col=c("#af8dc3","#7fbf7b");
LDAticks=c(-3, 0, 3);
lefse2$taxa=forcats::fct_rev(factor(lefse2$taxa));

#plot diverging plot
lefsep<-ggplot(lefse2, aes(x=taxa, 
                           y=LDA)) +
  geom_point(stat='identity', 
             aes(colour=region),
             size = 2)+
  scale_y_continuous(breaks=LDAticks)+
  scale_colour_manual(values=region_col,
                      breaks=c("Laikipia","Masai_Mara"))+
  facet_grid(~species, scales="free_x")+ 
  coord_flip() +
  theme_bw() + 
  labs(y="LDA Effect Size",
       x="", 
       colour="Enriched in") +
  theme(legend.position="bottom",
        legend.text=element_text(size=9),
        legend.title=element_text(size=9, face="bold"),
        axis.ticks = element_blank(),
        axis.text.y=element_text(size=6.5),
        axis.title.x=element_text(size=8, face="bold"),
        axis.text.x=element_text(size=8),
        strip.text = element_text(size =7.5, face="bold"))+
  guides(shape = guide_legend(override.aes = list(size = 1)));

plot(lefsep);

################################################################################
#             5. save Lefse plot
################################################################################
ggsave(filename="03_laikipia_mara_lefse_plot.pdf",
       device="pdf",path="./figures/Laikipia",
       plot=lefsep,
       width=7,
       height=4.5,
       units="in",
       dpi=500);