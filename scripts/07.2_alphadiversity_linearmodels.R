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

##CODE FOR: 
#A) running linear models testing whether gut microbiota alpha-diversity
#         varies with host family, species, and dietary guild
#B) running post-hoc comparisons if linear model was significant

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load sample metadata, and table of alpha diversity metrics
################################################################################
load("data/01_sample_metadata_filtered.Rdata")
load("data/07_alpha_values.Rdata")
alphadf=inner_join(alpha[,c(1,3,5,10)], meta[,c(1,2,5,8,12)], by="Group")


################################################################################
#             2. Run linear models of gut microbiota alpha-diversity ~
#     host dietary guild + host family, while controlling for sample month
#                               ALL HERBIVORES
################################################################################

#there are three alpha-diversity metrics of interest, so will use a for loop
#we will also save the model to a list for retrieval later

mymetrics=c("","Chao1", "Shannon","PD")
hlist <- list()  

for(i in 2:4)
{
  print(paste("Linear mixed model, all herbivores, for:", mymetrics[i]));
  hmod=lmer(alphadf[,i]~
              alphadf$diet_guild+
              alphadf$Family+
              (1|alphadf$sample_date),
            data=alphadf, na.action=na.fail);
  print(summary(hmod));
  print("Are host predictors statistically significant?");
  print(Anova(hmod));
  hlist[[i]]=hmod;
};

################################################################################
#             3. Run linear models of gut microbiota alpha-diversity ~
#     host dietary guild + host species, while controlling for sample month
#                               BOVIDS ONLY
################################################################################

#run the same model but on bovids only
bov=alphadf[alphadf$Family=="Bovidae",];
bov$species_short=factor(bov$species_short, 
                         levels=c("Buffalo","Cattle","Topi","Eland",
                                  "Impala","Gazelle","Dikdik"))

blist <- list()  
for(i in 2:4)
{
  print(paste("Linear mixed model, bovids only, for:", mymetrics[i]));
  bmod=lmer(bov[,i]~
              bov$diet_guild+
              bov$species_short+
              (1|bov$sample_date),
            data=bov, na.action=na.fail);
  print(summary(bmod));
  print("Are host predictors statistically significant?");
  print(Anova(bmod));
  blist[[i]]=bmod;
};


################################################################################
#             4. Run Tukey post-hoc comparisons on linear model: ALL HERBIVORES
################################################################################
#use for loop to run multiple comparison Tukey testing

#which diet groups differ from each other in terms of their gut microbiota 
#alphadiversity?
for(i in 2:4)
{
  print(paste("Tukey PostHoc Diet Comparisons, all herbivores, for:", 
              mymetrics[i]));
  print(summary(glht(hlist[[i]], 
                     linfct=mcp("alphadf$diet_guild"="Tukey")),
                test=adjusted("BH")));
};

#which host Families differ from each other in terms of their gut microbiota 
#alphadiversity?
for(i in 2:4)
{
  print(paste("Tukey PostHoc Host Family Comparisons, all herbivores, for:", 
              mymetrics[i]));
  print(summary(glht(hlist[[i]], 
                     linfct=mcp("alphadf$Family"="Tukey")),
                test=adjusted("BH")));
};


################################################################################
#             5. Run Tukey post-hoc comparisons on linear model: BOVIDS ONLY
################################################################################
#use for loop to run multiple comparison Tukey testing

#which diet groups differ from each other in terms of their gut microbiota 
#alphadiversity?
for(i in 2:4)
{
  print(paste("Tukey PostHoc Diet Comparisons, bovids only, for:", 
              mymetrics[i]));
  print(summary(glht(blist[[i]], 
                     linfct=mcp("bov$diet_guild"="Tukey")),
                test=adjusted("BH")));
};

#which host species differ from each other in terms of their gut microbiota 
#alphadiversity?
for(i in 2:4)
{
  print(paste("Tukey PostHoc Host Family Comparisons, all herbivores, for:", 
              mymetrics[i]));
  print(summary(glht(blist[[i]], 
                     linfct=mcp("bov$species_short"="Tukey")),
                test=adjusted("BH")));
};

##SPECIAL NOTE: If getting the error message "number of items to replace is not a multiple 
#of replacement length", for the code immediately above, rerun the bovid lmer models but put the
# host species term before the host diet term
