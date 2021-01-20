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

##CODE FOR: generating
#A) phylosymbiosis figure of gut microbiota similarity vs. host divergence time 

#B) phylosymbiosis dendogram based on gut microbiota similarity using hr clustering 

source(file="scripts/00_background.R"); #load necessary packages and specifications

################################################################################
#             1. Load matrix of host divergence times and matrix of
#                   microbiota distances
################################################################################

load("data/06_divergence_times_matrix.Rdata");
load("data/06_distances_mantel.Rdata");


################################################################################
#             2. Reorganize data frame to facilitate plotting in ggplot2
#                           ALL HERBIVORES
###############################################################################
#melt matrix so its structured as follows:
#host sp 1 | host sp 2 | divergence time | microbiota distance 
#code uses Bray-Curtis distances but can work for any distance matrix

#melt divergence time matrix
div=tibble::rownames_to_column(as.data.frame(divergence_times), 
                                "taxon_1");
div2=gather(div, key="taxon_2", 
             value="div_time", -taxon_1);

#melt microbiota distance matrix
micro=as.matrix(mantelbrayall);   #OR manteljacall, mantelunifracall, mantelunwunifracall
micro=tibble::rownames_to_column(as.data.frame(micro), 
                                 "taxon_1"); 
micro2=gather(micro, key="taxon_2", value="dissim", -taxon_1);

#merge the two data frames
all_herb<- div2 %>% right_join(micro2, by=c("taxon_1","taxon_2"))

#remove self comparisons (e.g. cattle vs cattle)
all_herb=all_herb[all_herb$dissim!=0,];


################################################################################
#             3. Plot scatterplot of gut microbiota dissimilarity 
#                       vs. host divergence time -- ALL HERBIVORES
###############################################################################
hps=ggplot(all_herb,aes(x=div_time, y=dissim))+
  geom_point(color="#969696", size=2)+
  geom_smooth(method=lm, color="black", se=FALSE)+
  theme_classic()+labs(title = "All herbivores")+
  labs(y="Microbiota dissimilarity (0-1)",
       x="Host Divergence time (mya)")+
  scale_x_continuous(breaks=seq(8,91, by=20))+
  theme(legend.title=element_blank(),text = element_text(size=12),
        legend.position="right",
        plot.title = element_text(size=12, face="bold"), 
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"));

plot(hps);

################################################################################
#             4. Reorganize data frame to facilitate plotting in ggplot2
#                           BOVIDS ONLY
###############################################################################
#subset divergence time matrix to only bovids
divergence_times_bov=divergence_times[c(1:4,7:8,10),c(1:4,7:8,10)];

#melt divergence time matrix
bovdiv=tibble::rownames_to_column(as.data.frame(divergence_times_bov), 
                                "taxon_1");
bovdiv2=gather(bovdiv, key="taxon_2", 
             value="div_time", -taxon_1);

#melt microbiota distance matrix
bovmicro=as.matrix(mantelbraybov);   #mantelbraybov, manteljacbov, mantelunifracbov, mantelunwunifracbov
bovmicro=tibble::rownames_to_column(as.data.frame(bovmicro), 
                                 "taxon_1");
bovmicro2=gather(bovmicro, key="taxon_2", value="dissim", -taxon_1);

#merge the two data frames
bov_herb<- bovdiv2 %>% right_join(bovmicro2, by=c("taxon_1","taxon_2"))

#remove self comparisons (e.g. cattle vs cattle)
bov_herb=bov_herb[bov_herb$dissim!=0,];

################################################################################
#             5. Plot scatterplot of gut microbiota dissimilarity 
#                       vs. host divergence time -- BOVIDS ONLY
###############################################################################
bps=ggplot(bov_herb,aes(x=div_time, y=dissim))+
  geom_point(color="#969696", size=2)+
  theme_classic()+labs(title = "Bovids only")+
  labs(y="",
       x="Host Divergence time (mya)")+
  scale_x_continuous(breaks=seq(8,14, by=2))+
  scale_y_continuous(breaks=c(0.42,0.6,0.76))+
  theme(legend.title=element_blank(),text = element_text(size=12),
        legend.position="right",
        plot.title = element_text(size=12, face="bold"), 
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"));

plot(bps);

################################################################################
#             6. Build gut microbiota dendogram of all study species
###############################################################################
#code uses microbiota Unifrac distances which seem more appropriate when dealing with phylogenetic distances
#but code works for any distance matrix

#apply hiearchical clustering to microbiota distance matrix
hr=hclust(mantelunifracall, method="average");
hr=as.phylo(hr);

#view numbered tree nodes
plot(hr,font=1); nodelabels(bg="white")

#rotate dendogram slightly
rt.14 <- rotate((hr), 12);
rt.14b <- rotate(rt.14, 13);
rt.14c <- rotate(rt.14b, 16);
rt.14d <- rotate(rt.14c, 17);
rt.14e <- rotate(rt.14d, 19);
rt.14f <- rotate(rt.14e, 20);
rt.14g <- rotate(rt.14f, 20); 
rt.14h <- rotate(rt.14g, 18);

#plot dendogram
plot.phylo(hr, type="p", 
           edge.col=1, edge.width=1.5, font=1,
           direction="leftwards",
           label.offset=0.005,show.node.label=TRUE, 
           no.margin=F);

##plot slightly rotated dendogram
plot.phylo(rt.14h, type="p", 
           edge.col=1, edge.width=1.5, font=1,
           direction="leftwards",
           label.offset=0.005,show.node.label=TRUE,
           show.tip.label=T,
           no.margin=F);

###Rotating tree does not affect clustering relationships, in rotated tree we:
#rotated the elephant vs. warthog/zebra polytomy
#rotated the gazelle, topi, impala polytomy
#moved the elephant,zebra, warthog node to the bottom


################################################################################
#             7. Save your figures
###############################################################################
#save your scatterplots into 1 figure
symbio=arrangeGrob(hps, bps, nrow=1);
ggsave(filename="06_phylosymbiosis_scatterplot.pdf",
       device="pdf",path="./figures",
       plot=symbio,
       width=8,
       height=3.5,
       units="in",
       dpi=500);

#save dendogram 
#manually save by going to Plot window: export > Save as PDF> save in "figures" folder
#saved as "06_microbiota_dendogram.pdf"