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

# find 2x and key features
# its only the 0-5 and 6-15
twos <- compart %>% 
  filter(str_detect(ForestType, "twice")) %>% 
  distinct(ForestType)


# organise
two_times <- compart %>% 
  filter(str_detect(ForestType, "00-05")|str_detect(ForestType, "06-15")) %>% 
  select(CompartmentName, LoggingRotation,Sp.01:Sp.27) %>%
  group_by(CompartmentName, LoggingRotation) %>%
  select_if(~!all(is.na(.))) %>% 
  summarise_all(.funs = function(x) (sum(x, na.rm = TRUE))) %>%
  mutate(LoggingRotation = case_when(
    LoggingRotation == 1 ~ "once",
    LoggingRotation == 2 ~ "twice")) %>% 
  ungroup()

two_times

# basic stats on once vs. twice
SR_per_compartment <- two_times %>% select_if(is.numeric) %>% 
  apply(., 1, sum)

shannon_per_compartment <- two_times %>% select_if(is.numeric) %>% 
  apply(., 1, diversity)

simpson_per_compartment <- two_times %>% select_if(is.numeric) %>% 
  apply(., 1, function(x) diversity(x, 'simpson'))

data_per_compartment <- 
  data.frame(select(two_times, CompartmentName, LoggingRotation), 
             SR_per_compartment, 
             shannon_per_compartment,
             simpson_per_compartment)

# iNEXT by compartment  
iNEXT_per_compartment <- two_times %>% select_if(is.numeric)

iNEXT_out <- list()

for(i in 1:15){
  tmp <- as.numeric(iNEXT_per_compartment[i,])
  use <- tmp[tmp>0]
  iNEXT_out[[i]] <- iNEXT(use)
}


# create data frame of once vs. twice 
# for basic iNEXT analysis
df <- two_times %>% 
  gather(Species, Richness, -CompartmentName, - LoggingRotation) %>% 
  group_by(Species, LoggingRotation) %>% 
  summarise(SR = sum(Richness))

once <- df %>% filter(LoggingRotation == "once") %>% select(Species, SR) %>% 
  rename(`Log Once` = SR) %>% ungroup()
twice <- df %>% filter(LoggingRotation == "twice") %>% select(Species, SR) %>% 
  rename(`Log Twice` = SR) %>% ungroup()

two_times_input <- data.frame(select(once, `Log Once`), select(twice, `Log Twice`))
rownames(two_times_input) <- once$Species
two_times_input

model <- iNEXT(two_times_input)
plot(model)
