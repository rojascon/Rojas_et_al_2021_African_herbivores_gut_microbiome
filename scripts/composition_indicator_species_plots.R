#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for plotting side barplots of bacterial species enriched in particular host groups
#bacterial species identified using indicator species analysis

source(file="scripts/background.R"); #load necessary packages and specifications

##load data from indicator species analyses
#see script "get_indicator_bacterial_sp.R"
hf=read.table("data/indicator_taxa_host_families.txt", sep="\t");
hg=read.table("data/indicator_taxa_host_dietguilds.txt", sep="\t");

############################### CLEAN UP DATA FILE; VERY MESSY ###################
#trim data frames to only pertinent data
hf=as.data.frame(hf[12:44,]); 
hg=as.data.frame(hg[10:41,]);

hf<- as.data.frame(hf[-c(2,10,17,24),]); 
colnames(hf)="datacol"; 
hg=as.data.frame(hg[-c(2,9,22),]); 
colnames(hg)="datacol"; 

#remove asterisks and trailing spaces from data
hf$A=gsub(pattern = "\\*", replacement="", 
          x=hf$datacol);

hg$A=gsub(pattern = "\\*", replacement="", 
          x=hg$datacol);

hf$B=gsub("\\s+", ",", gsub("^\\s+|\\s+$", "",hf$A));
hg$B=gsub("\\s+", ",", gsub("^\\s+|\\s+$", "",hg$A));

#make new data frames with 3 columns: 
#bacterial taxa | indicator species statistic | p-value
hfsplit=as.data.frame(str_split_fixed(hf$B, ",",3));
hgsplit=as.data.frame(str_split_fixed(hg$B, ",",3));

markhf= grep(pattern = "Group", x=as.character(hfsplit$V1));
markhg= grep(pattern = "Group", x=as.character(hgsplit$V1));

hfsplit$hostfam=c(
  rep("Giraffidae", markhf[2]-markhf[1]),
  rep("Suidae", markhf[3]-markhf[2]),
  rep("Equidae", markhf[4]-markhf[3]),
  rep("Elephantidae", nrow(hfsplit)-markhf[4]+1));

hgsplit$hostguild=c(
  rep("grazer", markhg[2]-markhg[1]),
  rep("browser", markhg[3]-markhg[2]),
  rep("mixed_feeder", nrow(hgsplit)-markhg[3]+1));

#remove unwanted rows
hfam<- as.data.frame(hfsplit[-markhf,]);
hguild<- as.data.frame(hgsplit[-markhg,]);

#rename columns of new data frames
colnames(hfam)=c("bacterial_family","stat","pvalue","host_family");
colnames(hguild)=c("bacterial_family","stat","pvalue","host_dietguild");

#coerce statistic column into numberical
hfam$stat=as.numeric(as.character(hfam$stat));
hguild$stat=as.numeric(as.character(hguild$stat));

#coerce bacterial taxa column to character
hfam$bacterial_family=as.character(hfam$bacterial_family);
hguild$bacterial_family=as.character(hguild$bacterial_family);

############################### PLOT BACTERIA TAXA ENRICHED IN HOST HERBIVORE FAMILIES ###################
#order levels of factor
hfam$host_family=factor(hfam$host_family, 
                        levels=c("Giraffidae","Suidae","Equidae","Elephantidae"));

#shorten bacterial taxa names
hfam$bacterial_family[hfam$bacterial_family=="Thermomicrobiales_JG30-KF-CM45"] = "Thermomicrobiales";

#select color palette
bar_col=c("#ffd92f","#1f78b4","#f781bf","palegreen");

#plot sideways barplot
ibarfam=ggplot(hfam, aes(x=bacterial_family, y=stat, fill=host_family))+ 
  geom_bar(stat="identity", width=0.5)+
  coord_flip()+
  scale_x_discrete(limits = hfam$bacterial_family)+
  scale_y_continuous(breaks=c(0, 0.25, 0.50, 0.75,1))+
  scale_fill_manual(values = bar_col)+
  labs(x = "",
       y = "Indicator Species statistic")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key.size=unit(0.6, "cm"),
        legend.text=element_text(size=7.2),
        legend.title=element_blank(),
        legend.position ="bottom",
        legend.spacing.y = unit(0, "cm"),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10));

##save image 
ggsave(filename="indicator_hostfamily.pdf",
       device="pdf",path="./figures",
       plot=ibarfam,
       width=5.3,
       height=4,
       units="in",
       dpi=500);

############################### PLOT BACTERIA TAXA ENRICHED IN HOST HERBIVORE DIET GUILDS ###################
#order levels of factor
hguild$host_dietguild=factor(hguild$host_dietguild, 
                        levels=c("grazer","browser","mixed_feeder"));

#shorten bacterial taxa names
hguild$bacterial_family[hguild$bacterial_family=="Coriobacteriales_Incertae_Sedis"] = "Coriobacteriales_InSed";

#select color palette
bar_col=c("#b3e2cd","#fdcdac","#bebada");

#plot sideways barplot
ibarguild=ggplot(hguild, aes(x=bacterial_family, y=stat, fill=host_dietguild))+ 
  geom_bar(stat="identity", width=0.5)+
  coord_flip()+
  scale_x_discrete(limits = hguild$bacterial_family)+
  scale_y_continuous(breaks=c(0, 0.25, 0.50, 0.75,1))+
  scale_fill_manual(values = bar_col)+
  labs(x = "",
       y = "Indicator Species statistic")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key.size=unit(0.6, "cm"),
        legend.text=element_text(size=9),
        legend.title=element_blank(),
        legend.position ="bottom",
        legend.spacing.y = unit(0, "cm"),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10));

##save image 
ggsave(filename="indicator_hostdietguild.pdf",
       device="pdf",path="./figures",
       plot=ibarguild,
       width=5.3,
       height=4,
       units="in",
       dpi=500);
