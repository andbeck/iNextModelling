## Nada Chapters 2 iNext Modelling ----

## libraries ----
library(iNEXT)
library(MASS) # stepAIC multiple regression tools
library(heplots) # effect sizes +
library(mvabund) # 
library(vegan)

library(gridExtra) # multi-panel ggplots
library(ggfortify) # diagnostics for models
library(ggrepel) # nice labelling of points
library(GGally)
library(broom)
library(ggalt)
library(tidyverse)

## Step 0: data import ----

# forest compartment
compart <- read_csv('./Chapter2_Compartments/DataSources/All_DTC3.csv')
habitat <- read_csv('./Chapter2_Compartments/DataSources/hab_summary_rot1_27.04.2018.csv')

# Step 1: calculate column sums by forest category and compartment name to get iNEXT vals ----

# isolate the data on logged 1 time or unlogged
# get category and compartment name and species
# sort and summarise_all creates column sums
# age a new column that combines age and name (e.g. the three compartments within each age category)

comp_age <- compart %>% 
  filter(!str_detect(ForestType, "twice")) %>%
  select(ForestCategory, CompartmentName, Sp.01:Sp.27) %>%
  group_by(ForestCategory, CompartmentName) %>%
  select_if(~!all(is.na(.))) %>% 
  summarise_all(.funs = function(x) (sum(x, na.rm = TRUE))) %>%
  mutate(Age_Name = paste(ForestCategory, CompartmentName, sep = ':')) %>% 
  mutate(ActualAge = case_when(CompartmentName %in% c("UTT35", "UTT41") ~ 2,
                               CompartmentName %in% c("PB87", "PS52") ~5,
                               CompartmentName == "JU93" ~ 6,
                               CompartmentName == "B12" ~ 10,
                               CompartmentName == "JU31Y" ~ 11,
                               CompartmentName == "JI52" ~ 20,
                               CompartmentName == "JI4" ~ 23,
                               CompartmentName == "B17" ~ 24,
                               CompartmentName %in% c("B18", "B19") ~ 30,
                               CompartmentName == "JI64" ~ 31,
                               CompartmentName %in% c("B18", "BT7") ~ 32,
                               CompartmentName == "JU31N" ~ 41,
                               CompartmentName %in% c("JU100", "PS26", "UTT27") ~ 0))

# Step 2: create Age based list for iNEXT ----

# for counting # 18 compartments
age_counter<-length(comp_age$CompartmentName)

# collection zone setup
compartment_input<- list()

# isolate the species count data for each age-compartment
for (i in seq_len(age_counter)){
  tmp <- as.data.frame(comp_age[i,]) # get a single compartment
  tmp <- select(tmp, -ForestCategory, - CompartmentName, -Age_Name) # isolate the species data
  tmp <- tmp[tmp>0] # clean away any 0's
  compartment_input[[i]]<-tmp # assign to the list above
}

# add sensible names to the list pieces
names(compartment_input) <- unique(comp_age$CompartmentName)

# Get Rid of entries with SR < 1 (iNEXT does not work without >2 species)
compartment_use <- compartment_input %>% keep(function(x) length(x)>1)
length(compartment_use)

# we are losing 3 because of only 1 species being found there
names(compartment_input)[!names(compartment_input) %in% names(compartment_use)]

# keep these for later adding (they won't work with iNEXT but we need the species richness = 1)
compartment_oneTransect <- compartment_input %>% keep(function(x) length(x)==1)

# set up category based data for iNEXT
category_use <- comp_age %>% 
  group_by(ForestCategory) %>% 
  select(starts_with("Sp.")) %>% 
  summarise_all(.funs = function(x) (sum(x, na.rm = TRUE)))

fc <- unique(category_use$ForestCategory)
rn <- colnames(category_use)[-1]
  
category_use <- data.frame(t(as.data.frame(category_use))[-1,])
names(category_use) <- c("five","fifteen","twentyfive","forty1","not")
category_use <- apply(category_use, 2, function(x) as.numeric(x))
rownames(category_use) <- rn
category_use <- data.frame(category_use)


#  Step 3: iNEXT modelling happens here -----

# run iNEXT (it works on a list)
# ignore warnings.
compartment_mod <- iNEXT(compartment_use, datatype = 'abundance', nboot = 999)
category_mod <- iNEXT(category_use)

# What we really want to do is to compute the estimate at a fixed min coverage
# at the minimum coverage (level  = NULL)
# set level = 0.8 for 80% for example (it will be extrapolated at anything beyond min)
coverage_min_diversity <- estimateD(compartment_use, base = 'coverage', level = 0.9) %>% 
  rename(Site = site) 

coverage_min_diversity_cat <- estimateD(category_use, base = 'coverage', level = 0.9) %>% 
  rename(Site = site) %>% filter(order == 0)

# evenness ----
## Species richness (S) and Pielou's evenness (J):
# S <- specnumber(BCI) ## rowSums(BCI > 0) does the same...
# J <- H/log(S)

even <- function(x){
  H <- diversity(x)
  S <- sum(x>0)
  J <- H/log(S)
  return(J)
}

# the pieces
compartment_use$PB87
diversity(compartment_use$PB87)
sum(compartment_use$PB87>0)

# test on single compartment
even(compartment_use$PB87)

# run on all compartments
evens <- t(map_df(compartment_input, function(x) even(x))) %>% data.frame()

# craft data frame for plotting
names(evens) <- "evenness"
evens <- rownames_to_column(evens)
names(evens)[1] <- "Compartment"
evens <- arrange(evens, Compartment)
habs <- arrange(habitat, Compartment)
even_habs <- data.frame(evenness = evens$evenness, habs)

e1 <- ggplot(even_habs, aes(x = Age.of.forest, y = evenness, label = Compartment))+
  geom_point(size = 3)+
  geom_text_repel()

e2 <- ggplot(even_habs, aes(x = herb_mean_height, y = evenness, label = Compartment))+
  geom_point(size = 3)+
  geom_text_repel()

e3 <- ggplot(even_habs, aes(x = litter_mean_depth, y = evenness, label = Compartment))+
  geom_point(size = 3)+
  geom_text_repel()

e4 <- ggplot(even_habs, aes(x = canopy_mean_closure, y = evenness, label = Compartment))+
  geom_point(size = 3)+
  geom_text_repel()

e5 <- ggplot(even_habs, aes(x = water_mean_number, y = evenness, label = Compartment))+
  geom_point(size = 3)+
  geom_text_repel()

e6<- ggplot(even_habs, aes(x = AdjacentUnlogged, y = evenness, label = Compartment))+
  geom_point(size = 3)+
  geom_text_repel()

grid.arrange(e1,e2,e3,e4,e5,e6, ncol = 3)


# visualise the results
# feel free to create more than one version of this!
par(mfrow = c(1,2))
plot(compartment_mod, type = 1)
plot(category_mod, type = 1)

ggplot(coverage_min_diversity_cat, aes(x = factor(Site, 
                                                  levels = c("five","fifteen","twentyfive","forty1","not")), 
                                       y = qD, 
                                       ymin = qD.LCL, ymax = qD.UCL,
                                       group = 1))+
  geom_point()+geom_line()+
  geom_errorbar(width =0)+
  xlab("Time Since Last Logging") + ylab("Species Richness") +
  theme_bw()

# Step 4: Collect Asymptotic estimates for downstream analyses and add habitat data ----
diversity_data_observations <- compartment_mod$AsyEst %>% arrange(as.character(Site))
diversity_data_coverage <- coverage_min_diversity %>% arrange(as.character(Site))

cbind(diversity_data_coverage$qD, diversity_data_observations$Estimator)

# organise the habitat data ----
# some compartments are not represented in the iNEXT data (no species)
habitat2 <- arrange(habitat, Compartment) %>% rename(Site = Compartment)
habitat_use <- filter(habitat2, Site %in% levels(diversity_data_coverage$Site))
names(habitat_use)

# get the various parts of iNEXT and habitat data together.----
# We need the habitat data represented three times each to match the 
# diversity data, which has estimates of SR, Shannon and Simpson for each compartment

df_prep <- data.frame(diversity_data_coverage, 
                 compartment_age = rep(habitat_use$Age.of.forest, each= 3),
                 canopy_closure = rep(habitat_use$canopy_mean_closure, each = 3),
                 leaf_litter_depth = rep(habitat_use$litter_mean_depth, each = 3),
                 understory_height = rep(habitat_use$herb_mean_height, each = 3),
                 No_WaterBodies = rep(habitat_use$water_mean_number, each = 3),
                 near_Primary = rep(habitat_use$AdjacentUnlogged, each = 3)) %>% 
  rename(coverage = SC) %>% 
  rename(diversity_estimate = qD) %>%
  select(-m,-method, -qD.LCL, -qD.UCL)

# ADD BACK IN THE Three Compartments with 1 SR
one_SR <- filter(habitat2, !Site %in% levels(diversity_data_coverage$Site))

one_DF <- bind_rows(replicate(3, one_SR, simplify = FALSE)) %>% 
  arrange(Site, Age.of.forest) %>% 
  mutate(diversity_estimate = 1, order = rep(0:2, 3), coverage = NA) %>% 
  select(Site, order, coverage, diversity_estimate, Age.of.forest, 
         contains('_mean_'), AdjacentUnlogged) %>% 
  as.data.frame()

names(one_DF) <- names(df_prep)

head(one_DF)
head(df_prep)

df <- rbind(df_prep, one_DF)
head(df)
# quick check
head(df)

## Exploratory Plotting of diversity and habitat data ----

# Isolate Species Richness data (order/q = 0)
# unlogged currenty as 200 yrs old.
SR_data <- filter(df, order == 0) %>% 
  mutate(compartment_age = as.numeric(as.character(compartment_age)))

Simp_data <- filter(df, order == 2) %>% 
  mutate(compartment_age = as.numeric(as.character(compartment_age)))

Shan_data <- filter(df, order == 1) %>% 
  mutate(compartment_age = as.numeric(as.character(compartment_age)))

head(SR_data)
tail(SR_data)

# pairs correlations
ggpairs(select(SR_data,compartment_age:near_Primary))

# fascinating groups....
SR_plot <- ggplot(SR_data, aes(x = compartment_age, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  labs(x = "Compartment Age (Years Since Last Logging)",
       y = "Species Richness Estimate") +
  theme_bw()

Simp_plot <- ggplot(Simp_data, aes(x = compartment_age, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  labs(x = "Compartment Age (Years Since Last Logging)",
       y = "Simpsons Diversity Estimate") +
  theme_bw()

grid.arrange(SR_plot, Simp_plot, ncol = 1)

# Formal plots ----
SR_age <- ggplot(SR_data, aes(x = compartment_age, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  labs(x = "Compartment Age (Years Since Last Logging)",
       y = "Species Richness Estimate") +
  theme_bw()

SR_canopy_closure <- ggplot(SR_data, aes(x = canopy_closure, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  labs(x = "Canopy Closure",
       y = "Species Richness Estimate") +
  theme_bw()

SR_leaf_litter_depth <- ggplot(SR_data, aes(x = leaf_litter_depth, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  labs(x = "Leaf Litter Depth",
       y = "Species Richness Estimate") +
  theme_bw()

SR_understory_height <- ggplot(SR_data, aes(x = understory_height, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  labs(x = "Understory Height",
       y = "Species Richness Estimate") +
  theme_bw()

SR_No_WaterBodies <- ggplot(SR_data, aes(x = No_WaterBodies, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  labs(x = "No.WaterBodies",
       y = "Species Richness Estimate") +
  theme_bw()

SR_Primary <- ggplot(SR_data, aes(x = near_Primary, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  labs(x = "near_Primary",
       y = "Species Richness Estimate") +
  theme_bw()

grid.arrange(SR_age, SR_canopy_closure, 
             SR_leaf_litter_depth, SR_understory_height, 
             SR_No_WaterBodies, SR_Primary,
             ncol = 3)

## linear modelling with iNEXT Asymptotes ----

# Species Richness
# linear model with all terms; polynomial 2 and 1 
mod_full_poly <- lm(diversity_estimate ~ poly(compartment_age,2) + 
                 poly(canopy_closure,2) + poly(leaf_litter_depth,2) + poly(understory_height,2) + 
                 No_WaterBodies + near_Primary, 
                 data = SR_data)

mod_full_lin <- lm(diversity_estimate ~ compartment_age + 
                 canopy_closure + leaf_litter_depth + understory_height + 
                 No_WaterBodies + near_Primary, data = SR_data)

# step 1 - no difference between models via SS/F
anova(mod_full_poly, mod_full_lin)

# step 2 - but the residuals of mod_full_lin are terrible
par(mfrow = c(2,2))
plot(mod_full_lin)

# step 3 - apply log() transformation to linear model
mod_full_lin_log <- lm(log(diversity_estimate) ~ compartment_age + 
                 canopy_closure + leaf_litter_depth + understory_height + 
                 No_WaterBodies + near_Primary, data = SR_data)

# step 4 diagnostics OK now
par(mfrow = c(2,2))
plot(mod_full_lin_log, add.smooth = FALSE)

# step 5: Anova - nothing is signficant
Anova(mod_full_lin_log)

# Step 6: Medium effect sizes....
eta <- etasq(mod_full_lin_log, anova = TRUE)


# Simpsons Diversity
# linear model with all terms; polynomial 2 and 1 
mod_full_poly_simp <- lm(diversity_estimate ~ poly(compartment_age,2) + 
                      poly(canopy_closure,2) + poly(leaf_litter_depth,2) + poly(understory_height,2) + 
                      No_WaterBodies + near_Primary, 
                    data = Simp_data)

mod_full_lin_simp <- lm(diversity_estimate ~ compartment_age + 
                     canopy_closure + leaf_litter_depth + understory_height + 
                     No_WaterBodies + near_Primary, data = Simp_data)

# step 1 - no difference between models via SS/F
anova(mod_full_poly_simp, mod_full_lin_simp)

# step 2: check diagnostics of simple model
par(mfrow = c(2,2))
plot(mod_full_lin_simp, add.smooth = FALSE)

# Step 3: Inference - the residials are fine without log for Simpsons
# and we find a singificant effect of compartment age!
# simpsons diversity increases with time since logging.
Anova(mod_full_lin_simp)
summary(mod_full_lin_simp)

# Step 4: effect sizes
# and time since logging has a large effect size too!
# and proximity to primary has a medium effect size.
eta_simp <- etasq(mod_full_lin_simp, anova = TRUE)
eta_simp


# Shannons Diversity
# linear model with all terms; polynomial 2 and 1 
mod_full_poly_shan <- lm(diversity_estimate ~ poly(compartment_age,2) + 
                           poly(canopy_closure,2) + poly(leaf_litter_depth,2) + poly(understory_height,2) + 
                           No_WaterBodies + near_Primary, 
                         data = Shan_data)

mod_full_lin_shan <- lm(diversity_estimate ~ compartment_age + 
                          canopy_closure + leaf_litter_depth + understory_height + 
                          No_WaterBodies + near_Primary, data = Shan_data)

# step 1 - no difference between models via SS/F
anova(mod_full_poly_shan, mod_full_lin_shan)

# step 2: check diagnostics of simple model
par(mfrow = c(2,2))
plot(mod_full_lin_shan, add.smooth = FALSE)

# Step 3: Inference - the residials are fine without log for Simpsons
# and we find a singificant effect of compartment age!
# simpsons diversity increases with time since logging.
Anova(mod_full_lin_shan)
summary(mod_full_lin_shan)

# Step 4: effect sizes
# and time since logging has a large effect size too!
# and proximity to primary has a medium effect size.
eta_shan <- etasq(mod_full_lin_shan, anova = TRUE)
eta_shan

effect_size_data <- data.frame(Term = rownames(eta), 
           Species_Richness = eta$`Partial eta^2`,
           Simpson = eta_simp$`Partial eta^2`,
           Shannon = eta_shan$`Partial eta^2`)

effect_size_plot <- na.omit(
  gather(effect_size_data, key = q, value = Partial_eta, -Term))

ggplot(effect_size_plot, aes(x = Term, y = Partial_eta, colour = q))+
  geom_lollipop(size = 1)+
  ylab(expression(paste("Partial ", eta^2)))+
  facet_wrap(~q) +
  geom_hline(yintercept = c(0.01,0.06,0.28), linetype = 'dashed')+
  theme_bw(base_size = 15)+
  guides(colour = FALSE)+
  coord_flip()

out_anovas <- na.omit(data.frame(
  q = rep(c("Species Richness", "Shannon", "Simpson"), each = 7),
  rbind(
    tidy(eta), tidy(eta_shan), tidy(eta_simp))))

out_anovas <- out_anovas %>% mutate_if(is.numeric, funs(round(., 3))) %>% 
  rename(Term = term, `F` = statistic, SS = sumsq)

write_csv(out_anovas, path = "./Chapter2_Compartments/DataSources/out_anovas.csv")

# MVABUND approach ----

# does community composition vary by habitat.

# prepare the abundance data
abund_dat <- comp_age %>% 
  ungroup() %>% 
  select(starts_with("Sp.")) %>% 
  mvabund()

# prepare the habitat data
hab_dat <- SR_data %>% 
  select(-Site, -order, -coverage, - diversity_estimate) %>% 
  as.data.frame()

# view some of the abundance data
plot(abund_dat, hab_dat$near_Primary)
plot(abund_dat, cut(hab_dat$compartment_age, 3))

# build the multivariate model of site x species abundance ~ habitat
mv_model <- manyglm(abund_dat ~ compartment_age+
                   canopy_closure +
                   leaf_litter_depth +
                   understory_height +
                   No_WaterBodies +
                   near_Primary, 
                 data = hab_dat, family = "negbinomial")

# diagnostics (some hint of mean-var with poisson.... negBin much better)
par(mfrow = c(1,3))
plot(mv_model, which = 1:3)

# model inference
set.seed(260518)
anova(mv_model)

# betapart to dig deeper?
library(betapart)
library(vegan)

beta_dat <- comp_age %>% 
  ungroup() %>% 
  select(starts_with("Sp."))

# create beta dissimilarity based on abundance
betaCore <- beta.pair.abund(beta_dat)
betaCore <- bray.part(beta_dat)
betaCore

# the MVabund model highlighted leaf litter and No.Water Bodies
# lets create a 3 level factor for each of these, and assess beta
timeSince_3fac <- cut(hab_dat$compartment_age, breaks = c(0,18,50,200), 
                     labels = c("recent", "recovered", "primary-like")) 
canopy_3fac <- cut(hab_dat$canopy_closure, 3, 
                   labels = c("closed", "speckled", "open"))
litter_3fac <- cut(hab_dat$leaf_litter_depth, 3, 
                   labels = c("shallow", "med", "deep"))
understory_3fac <- cut(hab_dat$understory_height, breaks = c(0,65,80,100),
                      labels = c("low", "med", "tall"))
water_3fac <- cut(hab_dat$No_WaterBodies, 3, 
                  labels = c("few", "mid", "many"))

par(mfrow = c(2,3))
plot(betadisper(betaCore[[3]], timeSince_3fac), main = "Time Since Logging")
plot(betadisper(betaCore[[3]], canopy_3fac), main = "Canopy Closure")
plot(betadisper(betaCore[[3]], litter_3fac), main = "Litter Depth **")
plot(betadisper(betaCore[[3]], understory_3fac), main = "Understory Height")
plot(betadisper(betaCore[[3]], water_3fac), main = "Water Bodies **")
plot(betadisper(betaCore[[3]], hab_dat$near_Primary), main = "Near Primary")

# betapart approach does not show any signficant changes
anova(betadisper(betaCore[[3]], timeSince_3fac))
anova(betadisper(betaCore[[3]], canopy_3fac))
anova(betadisper(betaCore[[3]], litter_3fac))
anova(betadisper(betaCore[[3]], understory_3fac))
anova(betadisper(betaCore[[3]], water_3fac))
anova(betadisper(betaCore[[3]], hab_dat$near_Primary))


# model averaging approach ----
library(MuMIn)

fullMod_dredge <- update(mod_full_lin, na.action = na.fail)
dd <- dredge(fullMod_dredge, subset = 
               dc(compartment_age, I(compartment_age^2)) && 
               dc(canopy_closure, I(canopy_closure^2)) &&
               dc(leaf_litter_depth, I(leaf_litter_depth^2)) && 
               dc(understory_height, I(understory_height^2)))
  

# set for delta AICc < 2 : the AICc is compared to the best model
# we want ones that are close to the best model
# delta = 'difference in'
summary(model.avg(dd, subset = delta < 3))

# best model
summary(get.models(dd, 1)[[1]])


