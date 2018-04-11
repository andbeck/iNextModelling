# elevation AccPres

# no compartment name or transect number

# ---------------------------------------
# setup with libraries

library(tidyverse)
library(iNEXT)

# ---------------------------------------
# source the function
# this grabs the function and loads it into R's brain
# to use.

source('Nada_Function_MAD.R')

# ---------------------------------------

HillDipt <- read.csv("./CrossVal_Example/HillDipt.csv")

glimpse(HillDipt)
unique(HillDipt$MountainName)
unique(HillDipt$Elevation)

# MODIFY DATA INPUT TO MATCH ACC PRES EXPECTATIONS
useDat <- HillDipt %>% 
  rename(CompartmentName = MountainName,TransectNo = Elevation) %>% 
  select(CompartmentName, TransectNo, Sp.01:Sp.28)

glimpse(useDat)
names(useDat)

# # This makes sure the start of AccPres() works
# wrk <- unite(useDat, col = Comp_Trans, CompartmentName,TransectNo, sep = ":")
# wrk <- wrk %>% select(Comp_Trans, starts_with("Sp."))
# 
# Trans<-wrk %>% 
#   # select(36,9:35) %>% # no vis count
#   group_by(Comp_Trans) %>%
#   summarise_all(funs(sum), na.rm = TRUE) # adds up all transects within compartment - date
# 
# AllData <- Trans %>% select(-Comp_Trans) %>% 
#   summarise_all(funs(sum)) %>% 
#   as.numeric()
# 
# Master <- iNEXT(AllData)
# plot(Master)

# ---------------------------------------
# Use the function on the  data frame to retrieve 
# a table of Precision and MAE estimates
# This will take ~ 30 seconds
# it will print some errors too... 

AccPres(useDat)

