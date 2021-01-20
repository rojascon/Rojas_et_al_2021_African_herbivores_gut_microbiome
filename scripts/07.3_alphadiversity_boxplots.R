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

##CODE FOR: generating side by side boxplots of gut microbiota alpha diversity
#just plotting Shannon alpha diversity, but can plot any alpha-diversity metric

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load sample metadata, and table of alpha diversity metrics
################################################################################
load("data/01_sample_metadata_filtered.Rdata")
load("data/07_alpha_values.Rdata")
alphadf=inner_join(alpha[,c(1,5)], meta[,c(1,5,8,12)], by="Group")

################################################################################
#             2. rearrange data frames to facilitate plotting in ggplot2
#                            ALL HERBIVORES
################################################################################
#melt alphadf for host family boxplots (A)
#melt alphadf for host dietary guild boxplots (B)
#then merge the two melted df(C)
#need to do this because essentially plotting the same data in both boxplots

A=alphadf[,c(1,2,4)]; 
A$cat="Host family"; 
colnames(A)=c("Group","value","variable","category") ;

B=alphadf[,c(1,2,5)]; 
B$cat="Host dietary guild";
colnames(B)=c("Group","value","variable","category");

C=rbind(A,B);
C$category=factor(C$cat, levels=c("Host family","Host dietary guild"));

################################################################################
#             3. plot boxplots of alpha diversity color coded by
#                   host family and host dietary guild
#                       ALL HERBIVORES
################################################################################
#select color palette
my_col=c("#fb6a4a","#ffd92f","#1f78b4","#f781bf","#74c476","#b3e2cd","#fdcdac","#bebada");

#create plot
herb_box=ggplot(data=C, 
                mapping=aes(x=variable,y=value, fill=variable))+
  geom_boxplot()+
  facet_wrap(~category, scales="free_x",nrow = 1)+ 
  theme_bw()+ 
  labs(x = "",
       y = "Shannon Diversity",
       title="All hervibores")+
  ylim(3.4,7)+
  scale_fill_manual(values=my_col)+
  theme(legend.position="none", 
        plot.title = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x=element_text(size=10, angle = 45, vjust = 0.66),
        strip.text = element_text(size =11, face="bold"));

plot(herb_box);

################################################################################
#             4. rearrange data frames to facilitate plotting in ggplot2
#                            BOVIDS ONLY
################################################################################

#subset alphaf to only include bovids
bov=alphadf[alphadf$Family=="Bovidae",]; 
bov$species_short=factor(bov$species_short, 
                         levels=c("Buffalo","Cattle","Topi","Eland",
                                  "Impala","Gazelle","Dikdik"))

#melt bov for host species boxplots (D)
#melt bov for host dietary guild boxplots (E)
#then merge the two melted df(G)
#need to do this because essentially plotting the same data in both boxplots

D=bov[,c(1,2,3)]; 
D$cat="Host species"; 
colnames(D)=c("Group","value","variable","category") ;

E=bov[,c(1,2,5)]; 
E$cat="Host dietary guild";
colnames(E)=c("Group","value","variable","category");

G=rbind(D,E);
G$category=factor(G$cat, levels=c("Host species","Host dietary guild"));


################################################################################
#             3. plot boxplots of alpha diversity color coded by
#                   host species and host dietary guild
#                       BOVIDS ONLY
################################################################################
#select color palette
my_col2=c("#1b9e77", "#d95f02", "darkorchid2", "#e7298a","steelblue1", 
          "#fbb4ae","#e6ab02","#b3e2cd","#fdcdac","#bebada");

#create plot
bov_box=ggplot(data=G, 
               mapping=aes(x=variable,y=value, fill=variable))+
  geom_boxplot()+
  facet_wrap(~category, scales="free_x",nrow = 1)+ 
  theme_bw()+ 
  labs(x = "",
       y = "Shannon Diversity",
       title="Bovids only")+
  ylim(3.4,7)+
  scale_fill_manual(values=my_col2)+
  theme(legend.position="none", 
        plot.title = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x=element_text(size=10, angle = 45, vjust=0.66),
        strip.text = element_text(size =11, face="bold"));

plot(bov_box);


################################################################################
#             6. Save your boxplots
################################################################################
##save all herbivore boxplot
ggsave(filename="07_alphadiv_boxplot_allherbivores.pdf",
       device="pdf",path="./figures",
       plot=herb_box,
       width=6.5,
       height=4,
       units="in",
       dpi=500);

##save bovids only boxplot
ggsave(filename="07_alphadiv_boxplot_bovids.pdf",
       device="pdf",path="./figures",
       plot=bov_box,
       width=6.5,
       height=4,
       units="in",
       dpi=500);