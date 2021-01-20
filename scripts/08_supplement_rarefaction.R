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

##CODE FOR: generating rarefaction curve of ASV richness after samples were
#subsampled to 17,000 sequences for alpha-diversity analyses

#NOTE: there will be a bit of a lag to run code from beginning to end
#because of large amount of data displayed

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Generate rarefaction curve data using mothur software                 
################################################################################
#Step1: open data file "data/07_output.mothur.txt" in Excel and save as .txt
#for some reason, mothur insists on this step ^^^

#Step2: open mothur and run command:
#       rarefaction.single(shared=07_output.mothur.txt, calc=sobs, freq=1)
#Code will take a few minutes to run

#Step3: rename output file as "08_mothur_rarefaction.txt" and place in the 
#data directory

################################################################################
#             2. Clean and arrange data for ggplot2                 
################################################################################
#read in data output from mothur after rarefaction
mrare=read.table("data/08_mothur_rarefaction.txt", 
                              sep="\t", header=TRUE);

#load sample metadata
load("data/01_sample_metadata_filtered.Rdata")

#clean up rarefaction curve data
raref=mrare[, grepl("X0.03.*",names(mrare))]; 
colnames(raref)=gsub("X0.03.", "", colnames(raref)); 

raref=cbind(mrare$numsampled,raref);
colnames(raref)[1]="numsampled";

#melt table for ggplot2 (timepoint|sample|numOTUs)
rf<- reshape2::melt(raref, id.vars="numsampled",value.name = "ASV_Richness");
colnames(rf)=c("TimeSampled", "Group", "ASV_Richness");
rare_meta=merge(rf, meta, by="Group");


################################################################################
#             3. Make plot of sample rarefaction curves                 
################################################################################

#set color-palette (1 color for each host species)
my_col=c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462",
         "#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f");

#set tick marks for x-axis 
reads=seq(0, nrow(mrare), by=2500);

#plot! 
rfc=ggplot(data=rare_meta) + 
  geom_line(mapping=aes(x=TimeSampled, 
                        y=ASV_Richness, 
                        col=species_short, 
                        group=Group))+ 
  scale_x_continuous(breaks=reads)+
  lims(y=c(0,17000))+
  scale_colour_manual(values=my_col)+
  labs(x = "Number of Reads",
       y = "ASV Richness",
       colour= "")+
  theme_bw() +
  guides(color=guide_legend(override.aes = list(size=3)))+
  theme(panel.background= element_rect(colour = "black", size=2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="right", 
        legend.text = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.x = element_text(size=14, face="bold"));

#plot(rfc);

################################################################################
#             4. save plot              
################################################################################
#NOTE: there is a bit of a lag because of large amount of data displayed
ggsave(filename="08_rarefaction_curves.pdf",
       device="pdf",path="./figures",
       plot=rfc,
       width=9,
       height=6,
       units="in",
       dpi=500);
