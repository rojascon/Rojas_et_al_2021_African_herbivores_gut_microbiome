#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for generating a table that shows the amount of samples collecter per herbivore species
#for each month (March-June)

source(file="scripts/background.R"); #load necessary packages and specifications

#load metadata with baboon samples removed, as well as samples that did not amplify
load("data/sample_metadata_beta.Rdata");

#make a dataframe of sample counts for each species for each month
dens=data.frame("March" = c(0, 13, 6, 2, 8, 3,29, 16,1,1,1), 
                "April" = c(9,0,1,1,3,1,2,1,0,0,5), 
                "May" = c(7,0,10,0,8,8,0,1,1,0,4),
                "June"=c(1,0,2,5,1,2,0,0,3,7,2))
dens$sp=c("Buffalo","Cattle","Topi","Eland","Impala","Gazelle",
                 "Dikdik","Giraffe","Zebra","Warthog","Elephant")

mbar<-reshape2::melt(dens,id.vars="sp", value.name = "count");
colnames(mbar)[2]="month"
mbar$sp=factor(mbar$sp, 
                levels=c("Buffalo","Cattle","Topi","Eland","Impala","Gazelle",
                 "Dikdik","Giraffe","Zebra","Warthog","Elephant"))

#save color gradient (1 color for each host species)
my_col=c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462",
         "#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f");

#create stacked
mbarp=ggplot(data=mbar, 
              mapping=aes(x=month,y=count, fill=sp))+
  geom_bar(stat ="identity")+
  theme_bw()+ 
  labs(x = "",
       y = "Number of Samples",
       fill="Host species")+
  scale_fill_manual(values=my_col)+
  #scale_fill_brewer(palette = "Paired")+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=13, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size=13, face="bold"),
        axis.title.x = element_text(size=13, face="bold"),
        axis.text.y = element_text(size=12),
        axis.text.x=element_text(size=12));

plot(mbarp);

##save image 
ggsave(filename="stacked_samples_bymonth.pdf",
       device="pdf",path="./figures",
       plot=mbarp,
       width=5.5,
       height=4,
       units="in",
       dpi=500);
