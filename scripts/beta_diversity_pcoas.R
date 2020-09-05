#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for constructing PCoAs based on beta-diversity distance matrices
#Bray-Curtis (bray), Jaccard (jac), Weighted Unifrac (wuni), Unweighted Unifrac (unwuni)

source(file="scripts/background.R"); #load necessary packages and specifications

#load metadata with baboon samples removed, as well as samples that did not amplify
load("data/sample_metadata_beta.Rdata");

#make a bovid only metadata file 
metabovdf=metadf[metadf$Family=="Bovidae",];

#load distance matrices computed from ASV abundances
#see script "get_beta_diversity_distances.R"
#one set of distance matrices for all herbivores; another set for bovids only; 
#a third set for grazers, browsers, mixed feeders
load("data/betadiv_dist_objects.Rdata")

#######################  ALL HERBIVORES  ##############################
#PCOA COLOR CODED BY HOST FAMILY AND HOST DIETARY GUILD
#make PCoA coordinates
pcoa_dec=cmdscale(bray.dist, eig=TRUE);  ##bray.dist, jac.dist, wuni.dist, unwuni.dist
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,metadf,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

#color-palette
fam_col=c("#d53e4f","#ffd92f","#1f78b4","#f781bf","palegreen"); 
guild_col=c("#b3e2cd","#fdcdac","#bebada");

#plot the PCoA-color coded by host family
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

#combine the two plots into one plot
pcoas_all=arrangeGrob(fam_pcoa,guild_pcoa, nrow=1);
plot(pcoas_all);

##save image 
ggsave(filename="betadiv_pcoas_allherbivores.pdf",
       device="pdf",path="./figures",
       plot=pcoas_all,
       width=12,
       height=4,
       units="in",
       dpi=500);

#######################  BOVIDS ONLY  ##############################
#PCOA COLOR CODED BY HOST SPECIES AND HOST DIETARY GUILD
#make PCoA coordinates
pcoa_dec=cmdscale(braybov.dist, eig=TRUE);  ##braybov.dist, jacbov.dist, wunibov.dist, unwunibov.dist
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,metabovdf,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

#color-palette
species_col=c("#1b9e77", "#d95f02", "darkorchid2", "#e7298a",
              "#1f78b4","#fbb4ae","#e6ab02","#b3e2cd") 
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

#combine the two plots into one plot
pcoas_bov=arrangeGrob(species_pcoa,guild_pcoa, nrow=1);
plot(pcoas_bov);

##save image 
ggsave(filename="betadiv_pcoas_bovids.pdf",
       device="pdf",path="./figures",
       plot=pcoas_bov,
       width=12,
       height=4,
       units="in",
       dpi=500);

#######################  DIETARY GUILDS PCoAs  ##############################
#ONE PCOA FOR EACH DIETARY GUILD; PCOA COLOR CODED BY HOST SPECIES

############## GRAZERS ONLY
#make PCoA coordinates
pcoa_dec=cmdscale(braygrazer.dist, eig=TRUE);  ##braygrazer.dist, jacgrazer.dist, wunigrazer.dist, unwunigrazer.dist
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,metadf,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

#color-palette
grazer_col=c("#a6cee3","#1f78b4","#d9f0a3","#238b45")

#plot the PCoA-color coded by host family
grazer_pcoa=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=species_short),
             size = 4,
             shape=21)+
  scale_fill_manual(values=grazer_col)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="Grazers")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="bottom",
        legend.text=element_text(size=13),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE));

############## BROWSERS ONLY
#make PCoA coordinates
pcoa_dec=cmdscale(braybrowser.dist, eig=TRUE);  ##braybrowser.dist, jacbrowser.dist, wunibrowser.dist, unwunibrowser.dist
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,metadf,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

#color-palette
browser_col=c("#fee08b","#a65628") 

#plot the PCoA-color coded by host family
browser_pcoa=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=species_short),
             size = 4,
             shape=21)+
  scale_fill_manual(values=browser_col)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="Browsers")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="bottom",
        legend.text=element_text(size=13),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE));

############## MIXED-FEEDERS ONLY
#make PCoA coordinates
pcoa_dec=cmdscale(braymixedfeed.dist, eig=TRUE);  ##braymixedfeed.dist, jacmixedfeed.dist, wunimixedfeed.dist, unwunimixedfeed.dist
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,metadf,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

#color-palette
mixedfeed_col=c("#1b9e77", "#d95f02", "darkorchid2","#1f78b4","#fbb4ae")

#plot the PCoA-color coded by host family
mixedfeed_pcoa=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=species_short),
             size = 4,
             shape=21)+
  scale_fill_manual(values=mixedfeed_col)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="Mixed-feeders")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="bottom",
        legend.text=element_text(size=13),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE));

#combine the three plots into one plot
pcoas_guilds=arrangeGrob(grazer_pcoa,browser_pcoa, mixedfeed_pcoa, nrow=1);
plot(pcoas_guilds);

##save image 
ggsave(filename="betadiv_pcoas_guilds.pdf",
       device="pdf",path="./figures",
       plot=pcoas_guilds,
       width=12,
       height=4.5,
       units="in",
       dpi=500);
