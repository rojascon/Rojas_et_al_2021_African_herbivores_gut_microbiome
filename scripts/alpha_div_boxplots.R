#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for generating side by side boxplots of gut microbiota alpha diversity

source(file="scripts/background.R"); #load necessary packages and specifications
#^^^also loads meta data

##load table of alpha diversity values generated using phyloseq and picante
#see script "get_alphadiversity_values.R" for how to generate
load("data/alphadiv_values.Rdata");

#merge alpha diversity values with sample metadata
alphameta=merge(alphadiv, meta, by="Group");

#remove baboon samples because not part of study
alphameta=alphameta[alphameta$species_short!="Baboon",]; 

#######################  ALL HERBIVORES  ##############################
#prepare data for side by side boxplots of 
#alpha diversity metric ~ host family
#alpha diversity metric ~ diet guild
side=alphameta[,which(colnames(alphameta) %in% 
                       c("Group","Shannon"))];

md=meta[,which(colnames (meta) %in% 
               c("Group","Family","diet_guild"))];

sidemd=merge(side, md, by="Group"); 

#create one data frame for host family boxplots (A)
#create one data frame for host dietary guild boxplots (B)
#then merge (C)

A=sidemd[,1:3]; 
A$cat="Host family"; 
colnames(A)=c("Group","value","variable","category") ;

B=sidemd[,c(1:2,4)]; 
B$cat="Host dietary guild";
colnames(B)=c("Group","value","variable","category");

C=rbind(A,B);
C$category=factor(C$cat, levels=c("Host family","Host dietary guild"));

#select color palette
#5 colors for family; 3 colors for dietary guild
my_col=c("#fb6a4a","#ffd92f","#1f78b4","#f781bf","#74c476","#b3e2cd","#fdcdac","#bebada");

#create tick marks of y-axis
alpha_breaks=seq(0, round(max(C$value)), by=0.5);

#create plot
herb_box=ggplot(data=C, 
                mapping=aes(x=variable,y=value, fill=variable))+
  geom_boxplot()+
  facet_wrap(~category, scales="free_x",nrow = 1)+ 
  theme_bw()+ 
  labs(x = "",
       y = "Shannon Diversity",
       title="All hervibores")+
  scale_y_continuous(breaks=alpha_breaks)+
  scale_fill_manual(values=my_col)+
  theme(legend.position="none", 
        plot.title = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x=element_text(size=10, angle = 45, vjust = 0.66),
        strip.text = element_text(size =11, face="bold"));

plot(herb_box);

##save image 
ggsave(filename="alphadiv_boxplot_allherbivores.pdf",
       device="pdf",path="./figures",
       plot=herb_box,
       width=6.5,
       height=4,
       units="in",
       dpi=500);

#######################  BOVIDS ONLY  ##############################
#subset data
bov=alphameta[alphameta$Family=="Bovidae",]; 

#prepare data for side by side boxplots of 
#alpha diversity metric ~ host family
#alpha diversity metric ~ diet guild

side2=bov[,which(colnames(bov) %in% 
                        c("Group","Shannon"))];

md2=meta[,which(colnames (meta) %in% 
                 c("Group","species_short","diet_guild"))];

sidemd2=merge(side2, md2, by="Group"); 

#create one data frame for host family boxplots (A2)
#create one data frame for host dietary guild boxplots (B2)
#then merge (C2)

A2=sidemd2[,1:3]; 
A2$cat="Host species"; 
colnames(A2)=c("Group","value","variable","category");

B2=sidemd2[,c(1:2,4)]; 
B2$cat="Host dietary guild";
colnames(B2)=c("Group","value","variable","category");

C2=rbind(A2,B2);
C2$category=factor(C2$cat, levels=c("Host species","Host dietary guild"));

#select color palette
#7 colors for family; 3 colors for dietary guild
my_col2=c("#1b9e77", "#d95f02", "darkorchid2", "#e7298a","#1f78b4",
          "#fbb4ae","#e6ab02","#b3e2cd","#fdcdac","#bebada");

#create tick marks of y-axis
alpha_breaks=seq(0, round(max(C$value)), by=0.5);

#create plot
bov_box=ggplot(data=C2, 
                mapping=aes(x=variable,y=value, fill=variable))+
  geom_boxplot()+
  facet_wrap(~category, scales="free_x",nrow = 1)+ 
  theme_bw()+ 
  labs(x = "",
       y = "Shannon Diversity",
       title="Bovids only")+
  scale_y_continuous(breaks=alpha_breaks)+
  scale_fill_manual(values=my_col2)+
  theme(legend.position="none", 
        plot.title = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x=element_text(size=10, angle = 45, vjust=0.66),
        strip.text = element_text(size =11, face="bold"));

plot(bov_box);

##save image 
ggsave(filename="alphadiv_boxplot_bovids.pdf",
       device="pdf",path="./figures",
       plot=bov_box,
       width=6.5,
       height=4,
       units="in",
       dpi=500);