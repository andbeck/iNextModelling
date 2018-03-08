# Use the function

# ---------------------------------------
# setup with libraries

library(tidyverse)
library(iNEXT)

# ---------------------------------------
# source the function
# this grabs the function and loads it into R's brain
# to use.

source('Nada_Function_MAE.R')

# ---------------------------------------
# get the data
# you can get all of them if you want
# just give them each a name
# IGNORE THIS: Five<-read.csv("CrossVal_Example/FiveYr.csv")

Five<-read.csv("FiveYr.csv")

# quick look at them
glimpse(Five)
unique(Five$CompartmentName)

# ---------------------------------------
# Use the function on the  data frame to retrieve 
# a table of Precision and MAE estimates
# This will take ~ 30 seconds
# it will print some errors too... 

AccPres(Five)
