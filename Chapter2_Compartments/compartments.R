## Nada Chapters 2 Once Logged Data Preparation ----
###  THIS SCRIPT IS NOT IN USE ANYMORE - deprecated to 
# ALL COMPARTMENT_iNEXT for basics
# AND to OneTime and TwoTimes Analysis

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

comp_age_df <- compart %>% 
  mutate(ForestAge = ifelse(is.na(ForestAge), 200, ForestAge)) %>% 
  filter(!str_detect(ForestType, "twice")) %>%
  select(ForestCategory, CompartmentName, ForestAge, Sp.01:Sp.27) %>% 
  # create data frame
  gather(Species, Count, -ForestCategory, - CompartmentName, - ForestAge) %>%
  # omit NA's
  na.omit() %>% 
  # convert Counts to numeric
  mutate(Count = as.numeric(Count)) %>% 
  # set up grouping
  group_by(ForestCategory, CompartmentName, ForestAge, Species) %>% 
  # get totals among transect by compartment
  summarise(Total = sum(Count, na.rm = TRUE)) %>% 
  # ungroup
  ungroup() %>% 
  # select the core data
  select(Species, CompartmentName, Total) %>% 
  # spread out for iNEXT
  # compartment columns, species rows
  spread(CompartmentName, Total) %>% 
  # replace NA's with 0 counts
  mutate_all(.funs = funs(ifelse(is.na(.), 0, .))) %>% 
  # convert to data frame
  as.data.frame()

comp_age_df

# which compartments did we loose for having no species
lost <- compart %>% 
  filter(!str_detect(ForestType, "twice")) %>%
  select(ForestCategory, CompartmentName, ForestAge, Sp.01:Sp.27) %>% 
  # create data frame
  gather(Species, Count, -ForestCategory, - CompartmentName, - ForestAge) %>% 
  # group at Compartment to see totals
  group_by(CompartmentName) %>% 
  mutate(Count = as.numeric(Count)) %>% 
  summarise(Total = sum(Count, na.rm = TRUE)) %>% 
  # which ones are 0.
  filter(Total == 0)

lost # B18 BT7 UTT41

# organise rows and select the compartment columns
rns <- comp_age_df$Species
rownames(comp_age_df) <- rns
comp_age_use <- select(comp_age_df, -Species)

glimpse(comp_age_use)

#  Step 3: iNEXT modelling happens here -----

# run iNEXT (it works on a list)
# ignore warnings.
compartment_mod <- iNEXT(comp_age_use, datatype = 'abundance', nboot = 999)
plot(compartment_mod)

# What we really want to do is to compute the estimate at a fixed min coverage
# at the minimum coverage (level  = NULL)
# set level = 0.8 for 80% for example (it will be extrapolated at anything beyond min)

coverage_min_diversity <- estimateD(comp_age_use, "abundance", 
                                    base = 'coverage', level = 0.8) %>% 
  rename(Site = site)

# Step 4: Collect Asymptotic estimates for downstream analyses and add habitat data ----
diversity_data_observations <- compartment_mod$AsyEst %>% 
  arrange(as.character(Site))
diversity_data_coverage <- coverage_min_diversity %>% 
  arrange(as.character(Site))

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
habitat_use$Age.of.forest <- ifelse(habitat_use$Age.of.forest==200, 
                                    unlogged_age, habitat_use$Age.of.forest)
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

