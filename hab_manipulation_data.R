# Create Nada Habitat Summary Data

library(tidyverse)
library(GGally)

# get the raw data
hab_master <- read_csv("./SpeciesAdultData/RawHabitatData.csv")

# make sure it makes sense
names(hab_master)
distinct(hab_master, Compartment)
distinct(hab_master, LoggingRotation)

# separate into Litter, Canopy and Herbaceous

hab_Litter <- hab_master %>% 
  gather(key = SampleLitter, value = LeafLitter,
         LeafLitter1:LeafLitter12) %>% 
  select(Compartment, `Age of forest`, LoggingRotation, 
         SampleLitter, LeafLitter)
head(hab_Litter)
distinct(hab_Litter, Compartment)

hab_Canopy <- hab_master %>% 
  gather(key = SampleCanopy, value = CanopyClosure,
         CanopyClosure1:CanopyClosure12) %>% 
  select(Compartment, `Age of forest`, LoggingRotation, 
         SampleCanopy, CanopyClosure)

hab_Herb <- hab_master %>% 
  gather(key = SampleHerb, value = HerbHeight,
       HerbHeight1:HerbHeight6) %>% 
  select(Compartment, `Age of forest`, LoggingRotation, 
         SampleHerb, HerbHeight)

dim(hab_Herb); dim(hab_Canopy); dim(hab_Litter)

# get summary stats and isolate Rotation = 1 and UNLOGGED

Herb_means <- hab_Herb %>% 
  group_by(Compartment, `Age of forest`, LoggingRotation) %>% 
  summarise(herb_mean_height = mean(HerbHeight, na.rm = TRUE),
            herb_sd_height = sd(HerbHeight, na.rm = TRUE),
            herb_median_height = median(HerbHeight, na.rm = TRUE),
            herb_n_height = sum(!is.na(HerbHeight))) %>% 
  filter(LoggingRotation == 1 | LoggingRotation == "unlogged") %>% 
  ungroup()
head(Herb_means)

Canopy_means <- hab_Canopy %>% 
  group_by(Compartment, `Age of forest`, LoggingRotation) %>% 
  summarise(canopy_mean_closure = mean(CanopyClosure, na.rm = TRUE),
            canopy_sd_closure = sd(CanopyClosure, na.rm = TRUE),
            canopy_median_closure = median(CanopyClosure, na.rm = TRUE),
            canopy_n_closure = sum(!is.na(CanopyClosure))) %>% 
  filter(LoggingRotation == 1 | LoggingRotation == 'unlogged') %>% 
  ungroup() %>% 
  select(canopy_mean_closure, canopy_sd_closure, 
         canopy_median_closure, canopy_n_closure)
head(Canopy_means)

Litter_means <- hab_Litter %>% 
  group_by(Compartment, `Age of forest`, LoggingRotation) %>% 
  summarise(litter_mean_depth = mean(LeafLitter, na.rm = TRUE),
            litter_sd_depth = sd(LeafLitter, na.rm = TRUE),
            litter_median_depth = median(LeafLitter, na.rm = TRUE),
            litter_n_depth = sum(!is.na(LeafLitter))) %>% 
  filter(LoggingRotation == 1 | LoggingRotation == "unlogged") %>% 
  ungroup() %>% 
  select(litter_mean_depth, litter_sd_depth, 
         litter_median_depth, litter_n_depth)
head(Litter_means)

hab_summary_rot1 <- data.frame(Herb_means, Litter_means, Canopy_means)
hab_summary_rot1 <- select(hab_summary_rot1,
                           Compartment,
                           Age.of.forest,
                           contains("_mean_"),
                           contains("_sd_"),
                           contains("_median_"),
                           contains("_n_")
                           )
head(hab_summary_rot1)

# some checkings
hab_summary_rot1 %>% 
  group_by(Compartment, Age.of.forest) %>% 
  summarise(n()) %>% 
  arrange(Age.of.forest)

# write_csv(hab_summary_rot1, path = "./SpeciesAdultData/hab_summary_rot1.csv")

pairDat <- select(hab_summary_rot1, contains("_mean_"))
ggpairs(pairDat)

mean_median <- select(hab_summary_rot1, contains("_mean_"), contains("_median_"))

