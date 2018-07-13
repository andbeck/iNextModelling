## Nada Chapters 2 Once Logged Data Preparation ----

## libraries ----
library(iNEXT)
library(MASS) # stepAIC multiple regression tools
library(heplots) # effect sizes +
library(mvabund) # 
library(vegan)

library(tidyverse)

library(GGally)
library(broom)
library(gridExtra) # multi-panel ggplots
library(ggfortify) # diagnostics for models
library(ggrepel) # nice labelling of points

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

# how many compartments have only 1 species? THREE
compartment_input %>% keep(function(x) length(x)==1)

# updated 2.015 iNEXT now produces estiamtes of qD = 1 when there is one species.
# none of the plotting works.
tmpMod <- iNEXT(compartment_input)

# Get Rid of entries with SR < 1 (estimateD is not working with n = 1 species)
# but we know the estiamte is 1.
compartment_use <- compartment_input %>% keep(function(x) length(x)>1)
length(compartment_use)

# # we are losing 3 because of only 1 species being found there
# names(compartment_input)[!names(compartment_input) %in% names(compartment_use)]

# # keep these for later adding (they won't work with iNEXT but we need the species richness = 1)
# compartment_oneTransect <- compartment_input %>% keep(function(x) length(x)==1)

# # set up category based data for iNEXT
# category_use <- comp_age %>% 
#   group_by(ForestCategory) %>% 
#   select(starts_with("Sp.")) %>% 
#   summarise_all(.funs = function(x) (sum(x, na.rm = TRUE)))
# 
# fc <- unique(category_use$ForestCategory)
# rn <- colnames(category_use)[-1]
#   
# category_use <- data.frame(t(as.data.frame(category_use))[-1,])
# names(category_use) <- c("five","fifteen","twentyfive","forty1","not")
# category_use <- apply(category_use, 2, function(x) as.numeric(x))
# rownames(category_use) <- rn
# category_use <- data.frame(category_use)


#  Step 3: iNEXT modelling happens here -----

# run iNEXT (it works on a list)
# ignore warnings.
compartment_mod <- iNEXT(compartment_use, datatype = 'abundance', nboot = 999)
plot(compartment_mod)
#category_mod <- iNEXT(category_use)

# extra things for Nada
#compartment_mod2 <- iNEXT(compartment_use, q = c(0,1,2), datatype = 'abundance', nboot = 999)
#ggiNEXT(compartment_mod2, type = 1, facet.var = "order")

# What we really want to do is to compute the estimate at a fixed min coverage
# at the minimum coverage (level  = NULL)
# set level = 0.8 for 80% for example (it will be extrapolated at anything beyond min)

coverage_min_diversity <- estimateD(compartment_use, "abundance", base = 'coverage', level = 0.9) %>% 
  rename(Site = site) 

coverage_min_diversity_cat <- estimateD(category_use, base = 'coverage', level = 0.9) %>% 
  rename(Site = site) %>% filter(order == 0)


# visualise the results
# feel free to create more than one version of this!
plot(compartment_mod, type = 1)

# # potentially interesting when classified
# ggplot(coverage_min_diversity_cat, aes(x = factor(Site, 
#                                                   levels = c("five","fifteen","twentyfive","forty1","not")), 
#                                        y = qD, 
#                                        ymin = qD.LCL, ymax = qD.UCL,
#                                        group = 1))+
#   geom_point()+geom_line()+
#   geom_errorbar(width =0)+
#   xlab("Time Since Last Logging") + ylab("Species Richness") +
#   theme_bw()

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

# adjust here the age of the unlogged: options are 200 or 100
unlogged_age <- 200
habitat_use$Age.of.forest <- ifelse(habitat_use$Age.of.forest==200, unlogged_age, habitat_use$Age.of.forest)
habitat_use$Age.of.forest

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

################ Here ends iNEXT and Habitat Data Processing for the 1x stuff.

