# Create Nada Habitat Summary Data

library(tidyverse)

hab_master <- read_csv("./SpeciesAdultData/RawHabitatData.csv")

names(hab_master)

hab_Litter <- hab_master %>% 
  gather(key = SampleLitter, value = LeafLitter,
         LeafLitter1:LeafLitter12) %>% 
  select(Compartment, `Age of forest`, LoggingRotation, SampleLitter, LeafLitter)
head(hab_Litter)

hab_Canopy <- hab_master %>% 
  gather(key = SampleCanopy, value = CanopyClosure,
         CanopyClosure1:CanopyClosure12) %>% 
  select(Compartment, `Age of forest`, LoggingRotation, SampleCanopy, CanopyClosure)

hab_Herb <- hab_master %>% 
  gather(key = SampleHerb, value = HerbHeight,
       HerbHeight1:HerbHeight6) %>% 
  select(Compartment, `Age of forest`, LoggingRotation, SampleHerb, HerbHeight)

all_hab <- data.frame(
  hab_Litter,
  select(hab_Canopy, SampleCanopy, CanopyClosure),
  select(hab_Herb, SampleHerb, HerbHeight)
)

dim(all_hab)    

