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
  mutate(CompartmentName = MountainName,TransectNo = Elevation)

# ---------------------------------------
# Use the function on the  data frame to retrieve 
# a table of Precision and MAE estimates
# This will take ~ 30 seconds
# it will print some errors too... 

AccPres(TwFivePlus)
