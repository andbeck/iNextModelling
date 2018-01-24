## Nada Chapters 2-3 iNext Modelling ----

## libraries ----
library(iNEXT)
library(tidyverse)

## data import ----
# elevation
elev <- read.csv('SpeciesAdultData/All-Elevation.csv')

# forest compartment
compart <- read.csv('SpeciesAdultData/All_DTC2.csv')

## Working with Elevation Data ----
names(elev)
head(elev)

# group elevations into bands ----
# this is Option 2 for grouping the data - Option 1 is by Forest Type
elev_banded <- mutate(elev,
                    band = 
                      ifelse(Elevation %in% 200:500, "a-200-500",
                             ifelse(Elevation %in% 550:1000, "b-550-1000",
                                    ifelse(Elevation %in% 1050:1500, "c-1050-1500", 'd-1550-1900'))))


## Create Column SUMS (by species) to prep for List for iNEXT ----

# colums by band to get iNEXT vals
elev_banded<-elev_banded %>% 
  select(-MountainName, -Elevation, -ForestType, -Date) %>%
  group_by(band) %>%
  summarise_all(.funs = function(x) (sum(x, na.rm = TRUE)))

#colsums by forest band to get iNEXT vals
elev_type<-elev %>% 
  select(-MountainName, -Elevation, -Date) %>%
  group_by(ForestType) %>%
  summarise_all(.funs = function(x) (sum(x, na.rm = TRUE)))

#colsums by mountain AND forest band to get iNEXT vals
elev_MF<-elev %>% 
  select(-Elevation, -Date) %>%
  group_by(MountainName, ForestType) %>%
  summarise_all(.funs = function(x) (sum(x, na.rm = TRUE))) %>%
  # create single column of mountain and type name
  mutate(MName_FType = paste(MountainName, ForestType, sep = ':'))



## Create Lists for iNEXT ----

# create BAND list for iNEXT
band_counter<-length(unique(elev_banded$band))

elev_BAND_input<- list()

for (i in 1:band_counter){
  tmp <- as.data.frame(elev_banded[i,-c(1)])
  tmp<-tmp[tmp>0]
  elev_BAND_input[[i]]<-tmp
}

names(elev_BAND_input)<-unique(elev_banded$band)

# create TYPE list for iNEXT
type_counter<-length(unique(elev_type$ForestType))

elev_TYPE_input<- list()

for (i in 1:type_counter){
  tmp <- as.data.frame(elev_type[i,-c(1)])
  tmp<-tmp[tmp>0]
  elev_TYPE_input[[i]]<-tmp
}

names(elev_TYPE_input)<-unique(elev_type$ForestType)

# create MOUNTAIN:TYPE list for iNEXT
# If we want to estimate lots 5 mountain x n bands and then
# do regressions....
MN_FT_counter<-length(unique(elev_MF$MName_FType))

elev_MNFT_input<- list()

for (i in 1:MN_FT_counter){
  tmp <- as.data.frame(elev_MF[i,-c(1,2,31)])
  tmp<-tmp[tmp>0]
  elev_MNFT_input[[i]]<-tmp
}

names(elev_MNFT_input)<-unique(elev_MF$MName_FType)

# Get Rid of entries with NO SR estimated
# elev_MNFT_input <- Filter(length, elev_MNFT_input)

# Get Rid of entries with SR < 1 (iNEXT does not work without >2 species)
elev_MNFT_input <- elev_MNFT_input[lapply(elev_MNFT_input, length)>1]


# Run and plot iNEXT ----
band_mod <- iNEXT(elev_BAND_input, datatype = 'abundance', nboot = 999)
type_mod <- iNEXT(elev_TYPE_input, datatype = 'abundance', nboot = 999)
MNFT_mod <- iNEXT(elev_MNFT_input, datatype = 'abundance', nboot = 999)

par(mfrow = c(1,3))
plot(band_mod, type = 1)
plot(type_mod, type = 1)
plot(MNFT_mod, type = 1)

## Working with Compartment Data ----

