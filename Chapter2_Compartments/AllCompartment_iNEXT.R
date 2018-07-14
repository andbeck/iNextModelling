library(tidyverse)
library(iNEXT)

# forest compartment
compart <- read_csv('./Chapter2_Compartments/DataSources/All_DTC3.csv')
glimpse(compart)
distinct(compart, CompartmentName) # 26 compartments.

# Step 1: calculate column sums by compartment into df
# to apply iNEXT

comp_age_df <- compart %>% 
  mutate(ForestAge = ifelse(is.na(ForestAge), 200, ForestAge)) %>% 
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

tbl_df(comp_age_df) # 23/26 compartments.

# organise rows and select the compartment columns
rns <- comp_age_df$Species
rownames(comp_age_df) <- rns
comp_age_use <- select(comp_age_df, -Species)

glimpse(comp_age_use)

# Step 2 - which have 0 species....
lost <- compart %>% 
  select(ForestCategory, CompartmentName, ForestAge, Sp.01:Sp.27) %>% 
  # create data frame
  gather(Species, Count, -ForestCategory, - CompartmentName, - ForestAge) %>% 
  # group at Compartment to see totals
  group_by(CompartmentName) %>% 
  mutate(Count = as.numeric(Count)) %>% 
  summarise(Total = sum(Count, na.rm = TRUE)) %>% 
  # which ones are 0.
  filter(Total == 0)

lost # B18 BT7 UTT41 have 0 species.

# Step 3 iNEXT

# Step A: iNEXT Asymp
asymps <- iNEXT(comp_age_use)$AsyEst
asymps
write_csv(asymps, path = "~/Desktop/asymps.csv")

# Step B: estimateDs
coverage90 <- estimateD(comp_age_use, "abundance", 
                        base = 'coverage', level = 0.90)
write_csv(asymps, path = "~/Desktop/coverage90.csv")
