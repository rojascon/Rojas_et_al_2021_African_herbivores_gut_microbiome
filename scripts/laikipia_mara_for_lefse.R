#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#MASAI MARA && LAIKIPIA DATA COMBINED

#Code for generating ASV tables for LEfSe analyses
#LEfSe will identify the ASVs that are differentially enriched in Masai Mara giraffes vs.
#Laikipia giraffes, for example
#Also code for plotting the output from LEfSe

source(file="scripts/background.R"); #load necessary packages and specifications

#load ASV table output by DADA2 && metadata
#samples that did not amplify & baboon samples were removed
#see script "laikipia_mara_get_ASV_table.R"
load("data/laikipia_mara_filtered_ASV_table.Rdata");
load("data/laikipia_mara_filtered_metadata.Rdata");

#load ASV taxonomy output by DADA2
tax=read.csv("data/laikipia_mara_ASV_taxonomy_silva_v132.csv", header=T); 

#give every ASV a name (e.g. ASV1|genus classification)
tax$label=paste("ASV",sep="",rownames(tax));
tax$asv_genus=paste(tax$label,sep="|",tax$Genus);
rownames(tax)=tax$ASV; 
tax$ASV=NULL;

#replace ASV labels with ASV_genus name
asvlef=asv.tbl;
colnames(asvlef)=tax$asv_genus[
  match(colnames(asvlef),rownames(tax))];

##transform ASV counts to ASV proportions
asvlef<- apply(asvlef, 1, function(i) (i/sum(i)));
colSums(asvlef) ;
nrow(metadf)==ncol(asvlef);

############### MAKE ASV ABUNDANCE TABLES FOR EACH HOST SPECIES ################
#samples are labelled by their host's geographic region
#only ASVs >0.01% average relative abundance for that species are included in data frame

myspecies=levels(metadf$species_short);

for(i in 1:length(myspecies))
{
  A=metadf$Group[metadf$species_short==myspecies[i]];
  B=as.data.frame(asvlef[, colnames(asvlef) %in% A]);
  colnames(B)=metadf$region[
    match(colnames(B),metadf$Group)];
  B$avg=rowMeans(B);
  B=B[B$avg>0.0001,];
  B$avg=NULL;
  assign(paste(myspecies[i]), B); 
};

#save output files from for-loop as .txt for LEfSe
mydfs=list(Buffalo, Cattle, Eland, Impala, Giraffe, Warthog, Zebra, Elephant);

for (i in 1:length(myspecies)){
  filename <- paste("data/LM_",myspecies[i],"_lefse.txt", sep="");
  write.table(mydfs[[i]], filename, col.names=NA,row.names=TRUE,sep="\t")
};

############### RUN LEFSE STATS ON THE GALAXY PLATFORM ################

#Go to LEfSe on Galaxy (https://huttenhower.sph.harvard.edu/galaxy/)
#Upload data frames by going to Get Data > upload file from your Computer

#Go to LEfSe > Format Data for LEfSe 
#Select whether the vectors (features and meta-data information) are listed in 
#rows or columns: ROWS
#Select which row to use as class: #1 region

#Alpha value for the factorial Kruskal-Wallis test among classes = 0.05
#rename all output files as "lefse_output_speciesname.txt"
#compile all results into one file for R = "LM_output_lefse_allspecies.csv"


################### PLOT LEFSE OUTPUT ########################
##this plot will visually show which ASVs are enriched in Masai Mara herbivores 
##or Laikipia herbivores

#read in LefSe output files
lefse=read.csv("data/LM_output_lefse_allspecies.csv", header=T);
lefse[,7:9]=str_split_fixed(lefse$ASV, "_",3); lefse$V9=NULL;

colnames(lefse)=c("ASV","mean","region","LDA","pvalue","species","label","taxa");

lefse$species=factor(lefse$species, levels=c("Buffalo", "Cattle", "Eland", 
                                                "Impala", "Giraffe", "Warthog", 
                                                "Zebra", "Elephant"));
#adjust pvalues for multiple comparisons
lefse$padjusted=p.adjust(lefse$pvalue, method="BH");
lefse=lefse[lefse$label!="NA",];

#see taxonomic distribution of enriched ASVs
table(lefse$taxa);

#only retain ASVs with LDA > 3.2 for plotting purposes
#you can set this cutoff as any value between 2 and 4.30
lefse2=lefse[lefse$padjusted<0.05,];
lefse2=lefse2[lefse2$LDA>3.2,];

lefse2$LDA[lefse2$region=="Laikipia"]=
  (lefse2$LDA[lefse2$region=="Laikipia"])*-1

#set up color and LDA tickmark labels for plot
region_col=c("#af8dc3","#7fbf7b");
LDAticks=c(-3, 0, 3);
lefse2$taxa=forcats::fct_rev(factor(lefse2$taxa));

#make the plot
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

#save plot
ggsave(filename="laikipia_mara_lefse_plot.pdf",
       device="pdf",path="./figures",
       plot=lefsep,
       width=7,
       height=4.5,
       units="in",
       dpi=500);
