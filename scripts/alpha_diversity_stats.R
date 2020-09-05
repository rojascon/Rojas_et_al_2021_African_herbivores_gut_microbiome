#Manuscript: Host phylogeny and host ecology structure the mammalian gut microbiota 
#at different taxonomic scales
#Rojas et al. 2020
#Code for generating statistics of gut microbiota alpha-diversity values

source(file="scripts/background.R"); #load necessary packages and specifications
#^^^also loads meta data

##load table of alpha diversity values generated using phyloseq and picante
#see script "get_alphadiversity_values.R" for how to generate
load("data/alphadiv_values.Rdata");

#merge alpha diversity values with sample metadata
alphameta=merge(alphadiv, meta, by="Group");

#remove baboon samples because not part of study
alphameta=alphameta[alphameta$species_short!="Baboon",]; 
#attach(alphameta);

#run linear models of alpha diversity metric~ host predictor variables
#Chao1 Richness, Shannon diversity, Phylogenetic diversity
all_herbivores=lmer(PD~  ##Chao1, Shannon, or PD
             diet_guild+
             Family+
             rain_two_weeks+
             tmax_two_weeks+
             tmin_two_weeks+
             (1|sample_month),
           data=alphameta, na.action=na.fail);

#run the same model but on bovids only
bov=alphameta[alphameta$Family=="Bovidae",];
bov$species_short=factor(bov$species_short, 
                         levels=c("Buffalo","Cattle","Topi","Eland",
                                  "Impala","Gazelle","Dikdik"))
  
bovids_only=lmer(PD~   ##Chao1, Shannon, or PD
                    diet_guild+
                    species_short+
                    rain_two_weeks+
                    tmax_two_weeks+
                    tmin_two_weeks+
                    (1|sample_month),
                  data=bov, na.action=na.fail);


#view results and rsquared value of linear models
print(summary(all_herbivores));
print(summary(bovids_only));

#test if predictors are statistically significant 
print("Are predictors significant? (all herbivores)");
print(Anova(all_herbivores));
print("Are predictors significant? (bovids only)");
print(Anova(bovids_only));

#conduct post-hoc comparisons with Bonferroni corrected pvalues
#all herbivores-model
print("Tukey PostHoc ~ Diet (all herbivores)");
print(summary(glht(all_herbivores, 
             linfct=mcp(diet_guild="Tukey")),
        test=adjusted("BH")));

print("Tukey PostHoc ~ host family (all herbivores)");
print(summary(glht(all_herbivores, 
                   linfct=mcp(Family="Tukey")),
              test=adjusted("BH")));

#bovids only model
print("Tukey PostHoc ~ Diet (bovids only)");
print(summary(glht(bovids_only, 
                   linfct=mcp(diet_guild="Tukey")),
              test=adjusted("BH")));

print("Tukey PostHoc ~ host species (bovids only)");
print(summary(glht(bovids_only, 
                   linfct=mcp(species_short="Tukey")),
              test=adjusted("BH")));


