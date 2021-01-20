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

##CODE FOR: constructing PCoAs based on beta-diversity distance matrices
#Bray-Curtis (bray), Jaccard (jac), Weighted Unifrac (wuni), Unweighted Unifrac (unwuni)

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load beta diversity distance matrices (e.g. dissimilarity matrices) 
#                         & filtered metadata table                 
################################################################################
load("data/01_sample_metadata_filtered.Rdata");
load("data/05_dissimilarity_distances_beta.Rdata")


################################################################################
#             2. make PCoA ordination -- ALL HERBIVORES  
#       points are color coded by host family or host dietary guild
#       uses Bray-Curtis distances, but any distance matrix works
################################################################################

#calculate coordinates for PCoA
pcoa_dec=cmdscale(bray.dist, eig=TRUE);  ##bray.dist, jac.dist, wuni.dist, unwuni.dist
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,meta,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

#set the color-palette, one to color code by family and another by diet guild
fam_col=c("#d53e4f","#ffd92f","#1f78b4","#f781bf","palegreen"); 
guild_col=c("#b3e2cd","#fdcdac","#bebada");

#plot the PCoA-color coded by host family first
fam_pcoa=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=Family),
             size = 4,
             shape=21)+
  scale_fill_manual(values=fam_col)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="All herbivores")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));

#plot the PCoA-color coded by host dietary guild
guild_pcoa=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=diet_guild),
             size = 4,
             shape=21)+
  scale_fill_manual(values=guild_col)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));

#place the two plots side by side
pcoas_all=arrangeGrob(fam_pcoa,guild_pcoa, nrow=1);
#plot(pcoas_all);

##save image 
ggsave(filename="05_betadiv_pcoas_allherbivores.pdf",
       device="pdf",path="./figures",
       plot=pcoas_all,
       width=12,
       height=4,
       units="in",
       dpi=500);


################################################################################
#             3. make PCoA ordination -- BOVIDS ONLY  
#       cpoints olor coded by host species or host dietary guild
#       uses Bray-Curtis distances, but any distance matrix works
################################################################################

#make bovid specific data-frame
metabov=meta[meta$Family=="Bovidae",]
metabov$species_short=factor(metabov$species_short, 
                             levels=c("Buffalo","Cattle","Topi","Eland","Impala",
                                      "Gazelle","Dikdik"));

#make PCoA coordinates
pcoa_dec=cmdscale(braybov.dist, eig=TRUE);  ##braybov.dist, jacbov.dist, wunibov.dist, unwunibov.dist
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,metabov,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

#set the color-palette, one to color code by family and another by diet guild
species_col=c("#1b9e77", "#d95f02", "darkorchid2", "#e7298a",
              "steelblue1","#fbb4ae","#e6ab02","#b3e2cd");
guild_col=c("#b3e2cd","#fdcdac","#bebada");

#plot the PCoA-color coded by host family
species_pcoa=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=species_short),
             size = 4,
             shape=21)+
  scale_fill_manual(values=species_col)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="Bovids only")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));

#plot the PCoA-color coded by host dietary guild
guild_pcoa=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=diet_guild),
             size = 4,
             shape=21)+
  scale_fill_manual(values=guild_col)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));

#combine the two plots into one plot
pcoas_bov=arrangeGrob(species_pcoa,guild_pcoa, nrow=1);
#plot(pcoas_bov);

##save image 
ggsave(filename="05_betadiv_pcoas_bovids.pdf",
       device="pdf",path="./figures",
       plot=pcoas_bov,
       width=12,
       height=4,
       units="in",
       dpi=500);

