# elevation AccPres

# note modifications to deal with
# no compartment name or transect number

# ---------------------------------------
# setup with libraries

library(tidyverse)
library(iNEXT)

# ---------------------------------------
# source the function
# this grabs the function and loads it into R's brain
# to use.  This has been updated to manage the 
# compartment/transect -> mountain/elevation issues
# and species with all NA's issues

source('Nada_Function_MAD.R')

## Import Data -----------------------------
# These are defined by forest types
impData <- read_csv("./CrossVal_Example/HillDipt.csv")

# check it through
glimpse(impData)
unique(impData$MountainName)
unique(impData$Elevation)

# MODIFY DATA INPUT TO MATCH ACC PRES EXPECTATIONS
# renaming mountain to compartment and elevation to transect
# grabbing only the columns needed
useDat <- impData %>% 
  rename(CompartmentName = MountainName,TransectNo = Elevation) %>% 
  select(CompartmentName, TransectNo, contains("Sp."))

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
