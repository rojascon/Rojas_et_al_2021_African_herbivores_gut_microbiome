#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for testing phylosymbiosis (mantel tests of host divergence time vs. 
#host microbiota dissimilarity) and generating accompanying scatterplots

source(file="scripts/background.R"); #load necessary packages and specifications

#load host species divergence times
#see script "get_host_divergence_times.R"
load("data/host_divergence_times.Rdata");

#load gut microbiota distance matrices built for this mantel test
#see script "get_mantel_distances.R"
#one set of distance matrices for all herbivores; another set for bovids only
load("data/mantel_dist_objects.Rdata");

######################## RUN MANTEL TESTS ##########################

#run mantel test---all herbivores
colnames(divergence_times)==
  colnames(as.matrix(mantelbrayall)); #mantelbrayall, manteljacall, mantelunifracall, mantelunwunifracall

mantel(divergence_times, mantelbrayall, 
       method="spear", 
       permutations=999); 


#run mantel test---bovids only
divergence_times_bov=divergence_times[c(1:4,7:8,10),c(1:4,7:8,10)];

colnames(divergence_times_bov)==
  colnames(as.matrix(mantelbraybov));

mantel(divergence_times_bov, mantelbraybov,
       method="spear", 
       permutations=999);

######################## PHYLOSYMBIOSIS SCATTERPLOTS ##########################
#############scatterplot -- all herbivores

#reorganize data frame so its structured as follows:
#host sp 1 | host sp 2 | divergence time | microbiota distance
#reflects pairwise comparisons

divt=tibble::rownames_to_column(as.data.frame(divergence_times), 
                                 "taxon_1");
divt2=gather(divt, key="taxon_2", 
          value="div_time", -taxon_1);

micro=as.matrix(mantelbrayall);   #mantelbrayall, manteljacall, mantelunifracall, mantelunwunifracall
micro=tibble::rownames_to_column(as.data.frame(micro), 
                                 "taxon_1"); 
micro2=gather(micro, key="taxon_2", value="dissimi", -taxon_1);

#check contents match before combining the two data frames
divt2$taxon_1==micro2$taxon_1; 
divt2$taxon_2==micro2$taxon_2;

allherb_phylo=divt2;
allherb_phylo$distance=micro2$dissimi; 

#remove self comparisons (e.g. cattle vs cattle)
allherb_phylo=allherb_phylo[allherb_phylo$distance!=0,];

#make phylosymbiosis scatterplot
hps=ggplot(allherb_phylo,aes(x=div_time, y=distance))+
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

######################## PHYLOSYMBIOSIS SCATTERPLOTS ##########################
#############scatterplot -- bovids only

#reorganize data frame so its structured as follows:
#host sp 1 | host sp 2 | divergence time | microbiota distance
#reflects pairwise comparisons

divt=tibble::rownames_to_column(as.data.frame(divergence_times_bov), 
                                "taxon_1");
divt2=gather(divt, key="taxon_2", 
             value="div_time", -taxon_1);

micro=as.matrix(mantelbraybov);   #mantelbraybov, manteljacbov, mantelunifracbov, mantelunwunifracbov
micro=tibble::rownames_to_column(as.data.frame(micro), 
                                 "taxon_1");
micro2=gather(micro, key="taxon_2", value="dissimi", -taxon_1);

#check contents match before combining the two data frames
divt2$taxon_1==micro2$taxon_1; 
divt2$taxon_2==micro2$taxon_2;

bov_phylo=divt2;
bov_phylo$distance=micro2$dissimi; 

#remove self comparisons (e.g. cattle vs cattle)
bov_phylo=bov_phylo[bov_phylo$distance!=0,];

#make phylosymbiosis scatterplot
bps=ggplot(bov_phylo,aes(x=div_time, y=distance))+
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

########SAVE YOUR FIGURE
symbio=arrangeGrob(hps, bps, nrow=1); 

ggsave(filename="phylosymbiosis_scatterplot.pdf",
       device="pdf",path="./figures",
       plot=symbio,
       width=7.5,
       height=3.5,
       units="in",
       dpi=500);

