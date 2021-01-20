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
#working from Laikipia subdirectories within the main directories

##CODE FOR: constructing PCoA ordinations based on beta-diversity distance matrices
#Bray-Curtis (bray), Jaccard (jac), Weighted Unifrac (wuni), Unweighted Unifrac (unwuni)

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load beta diversity distance matrices (e.g. dissimilarity matrices) 
#                         & filtered metadata table                 
################################################################################
load("data/Laikipia/01_LM_sample_metadata_filtered.Rdata");
load("data/Laikipia/04_LM_dissimilarity_distances_beta.Rdata")


################################################################################
#             2. generate coordinates for PCoA ordination
################################################################################

#make PCoA coordinates
pcoa_dec=cmdscale(bray.dist, eig=TRUE);  ##bray.dist, jac.dist, wuni.dist, unwuni.dist
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,meta,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);


################################################################################
#             3. Plot PCoA ordination color coded by host species and 
#                  using shape to indicate host geographic region
################################################################################

#color-palette
my_col=c("#1b9e77", "#d95f02", "#e7298a","#225ea8",
         "#ffd92f","#f781bf","#43a2ca","palegreen"); 

#plot the PCoA-color coded by host family, and shape indicates region
LM_pcoa=ggplot(pcoa_met, aes(Axis1,Axis2, 
                             shape=factor(region, 
                                          labels=c("Masai Mara","Laikipia"))))+
  geom_point(aes(colour=species_short), size = 2)+
  scale_colour_manual(values=my_col)+
  scale_shape_manual(values = c(19,21))+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       colour="",
       title="Masai Mara and Laikipia herbivores")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_blank(),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));

#plot(LM_pcoa);

################################################################################
#             4. save PCoA ordination
################################################################################
ggsave(filename="04_laikipia_mara_pcoa.pdf",
       device="pdf",path="./figures/Laikipia",
       plot=LM_pcoa,
       width=6,
       height=4,
       units="in",
       dpi=500);