## Nada Chapters 2-3 iNext Modelling ----

## libraries ----
library(iNEXT)
library(MASS) # stepAIC multiple regression tools
library(heplots) # effect sizes +

library(tidyverse)
library(gridExtra) # multi-panel ggplots
library(ggfortify) # diagnostics for models
library(ggrepel) # nice labelling of points
library(GGally)
## Step 0: data import ----

# forest compartment
compart <- read.csv('SpeciesAdultData/All_DTC2.csv')
habitat <- read.csv('SpeciesAdultData/DTC-AverageHabitatCompt-Age(12.04.18).csv')

str(compart)
str(habitat)

habitat %>% 


# Step 1: calculate column sums by forest category and compartment name to get iNEXT vals ----

# isolate the data on logged 1 time or unlogged
# get category and compartment name and species
# sort and summarise_all creates column sums
# age a new column that combines age and name (e.g. the three compartments within each age category)

comp_age <- compart %>% 
  filter(!str_detect(ForestType, "twice")) %>%
  select(ForestCategory, CompartmentName, Sp.01:Sp.27) %>%
  group_by(ForestCategory, CompartmentName) %>%
  summarise_all(.funs = function(x) (sum(x, na.rm = TRUE))) %>%
  mutate(Age_Name = paste(ForestCategory, CompartmentName, sep = ':')) %>% 
  mutate(ActualAge = case_when(CompartmentName %in% c("UTT35", "UTT41") ~ 2,
                               CompartmentName %in% c("PB87", "PS52") ~5,
                               CompartmentName == "JU93" ~ 6,
                               CompartmentName == "B12" ~ 10,
                               CompartmentName == "JU31Y" ~ 11,
                               CompartmentName == "JI52" ~ 20,
                               CompartmentName == "JI4" ~ 23,
                               CompartmentName == "B17" ~ 24,
                               CompartmentName %in% c("B18", "B19") ~ 30,
                               CompartmentName == "JI64" ~ 31,
                               CompartmentName %in% c("B18", "BT7") ~ 32,
                               CompartmentName == "JU31N" ~ 41,
                               CompartmentName %in% c("JU100", "PS26", "UTT27") ~ 0))

# Step 2: create Age based list for iNEXT ----

# for counting # 18 compartments
age_counter<-length(comp_age$CompartmentName)

# collection zone setup
compartment_input<- list()

# isolate the species count data for each age-compartment
for (i in seq_len(age_counter)){
  tmp <- as.data.frame(comp_age[i,]) # get a single compartment
  tmp <- select(tmp, -ForestCategory, - CompartmentName, -Age_Name) # isolate the species data
  tmp <- tmp[tmp>0] # clean away any 0's
  compartment_input[[i]]<-tmp # assign to the list above
}

# add sensible names to the list pieces
names(compartment_input)<-unique(comp_age$CompartmentName)

# Get Rid of entries with SR < 1 (iNEXT does not work without >2 species)
compartment_use <- compartment_input %>% keep(function(x) length(x)>1)
length(compartment_use)

#  Step 3: iNEXT modelling happens here -----

# run iNEXT (it works on a list)
# ignore warnings.
compartment_mod <- iNEXT(compartment_use, datatype = 'abundance', nboot = 999)

# What we really want to do is to compute the estimate at a fixed min coverage
# at the minimum coverage (level  = NULL)
# set level = 0.8 for 80% for example (it will be extrapolated at anything beyond min)
coverage_min_diversity <- estimateD(compartment_use, base = 'coverage', level = NULL) %>% 
  rename(Site = site) 

# visualise the results
# feel free to create more than one version of this!
par(mfrow = c(1,2))
plot(compartment_mod, type = 1)
plot(compartment_mod, type = 3)

# Step 4: Collect Asymptotic estimates for downstream analyses and add habitat data ----
diversity_data_observations <- compartment_mod$AsyEst %>% arrange(as.character(Site))
diversity_data_coverage <- coverage_min_diversity %>% arrange(as.character(Site))
  
# organise the habitat data ----
# some compartments are not represented in the iNEXT data (no species)
habitat2 <- arrange(habitat, Compartment) %>% rename(Site = Compartment)
habitat_use <- filter(habitat2, Site %in% levels(diversity_data_coverage$Site))
names(habitat_use)

# get the various parts of iNEXT and habitat data together.----
# We need the habitat data represented three times each to match the 
# diversity data, which has estimates of SR, Shannon and Simpson for each compartment

df_prep <- data.frame(diversity_data_coverage, 
                 compartment_age = rep(habitat_use$ForestAge, each= 3),
                 canopy_closure = rep(habitat_use$CanopyClosure, each = 3),
                 leaf_litter_depth = rep(habitat_use$LeafLitterDepth, each = 3),
                 understory_height = rep(habitat_use$HerbPlantHeight, each = 3),
                 No.WaterBodies = rep(habitat_use$No.WaterBodies, each = 3),
                 WaterDistance = rep(habitat_use$WaterDistance, each = 3)) %>% 
  rename(coverage = SC) %>% 
  rename(diversity_estimate = qD) %>%
  select(-m,-method, -qD.LCL, -qD.UCL)

# ADD BACK IN THE Three Compartments with 0 SR
zeroSR <- filter(habitat2, !Site %in% levels(diversity_data_coverage$Site))

ZeroDF <- bind_rows(replicate(3, zeroSR, simplify = FALSE)) %>% 
  arrange(Site, ForestAge) %>% 
  mutate(diversity_estimate = 0, order = rep(0:2, 4), coverage = NA) %>% 
  select(Site, order, coverage, diversity_estimate, ForestAge, CanopyClosure:WaterDistance)
names(ZeroDF) <- names(df_prep)

head(ZeroDF)
head(df_prep)

df <- bind_rows(df, ZeroDF)

# quick check
head(df)

## Exploratory Plotting of diversity and habitat data ----

# Isolate Species Richness data (order/q = 0)
# temporarily remove Unlogged?  (set at age 0 above)
SR_data <- filter(df, order == 0) %>% filter(compartment_age!='Unlogged') %>% 
  mutate(compartment_age = as.numeric(as.character(compartment_age)))

SR_unlogged <- filter(df, order == 0) %>% filter(compartment_age =='Unlogged')

Simp_data <- filter(df, order == 2) %>% filter(compartment_age!='Unlogged') %>% 
  mutate(compartment_age = as.numeric(as.character(compartment_age)))

Simp_unlogged <- filter(df, order == 2) %>% filter(compartment_age =='Unlogged')

head(SR_data)
tail(SR_data)

# pairs correlations
ggpairs(select(SR_data,compartment_age:No.WaterBodies))

# fascinating groups....
SR_plot <- ggplot(SR_data, aes(x = compartment_age, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  geom_point(data = SR_unlogged, aes(x = 50, y = diversity_estimate), colour = 'red', size = 5)+
  labs(x = "Compartment Age (Years Since Last Logging)",
       y = "Species Richness Estimate") +
  annotate('text', x = 48, y = 3, angle = 90, label = "Unlogged", col = 'red')+
  theme_bw()

Simp_plot <- ggplot(Simp_data, aes(x = compartment_age, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  geom_point(data = Simp_unlogged, aes(x = 50, y = diversity_estimate), colour = 'red', size = 5)+
  labs(x = "Compartment Age (Years Since Last Logging)",
       y = "Simpsons Diversity Estimate") +
  annotate('text', x = 48, y = 2, angle = 90, label = "Unlogged", col = 'red') +
  theme_bw()

grid.arrange(SR_plot, Simp_plot)

# habitat plots ----
SR_canopy_closure <- ggplot(SR_data, aes(x = canopy_closure, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  labs(x = "Canopy Closure",
       y = "Species Richness Estimate") +
  theme_bw()

SR_leaf_litter_depth <- ggplot(SR_data, aes(x = leaf_litter_depth, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  labs(x = "Leaf Litter Depth",
       y = "Species Richness Estimate") +
  theme_bw()

SR_understory_height <- ggplot(SR_data, aes(x = understory_height, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  labs(x = "Understory Height",
       y = "Species Richness Estimate") +
  theme_bw()

SR_No.WaterBodies <- ggplot(SR_data, aes(x = No.WaterBodies, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=Site), point.padding = 0.5) +
  labs(x = "No.WaterBodies",
       y = "Species Richness Estimate") +
  theme_bw()

grid.arrange(SR_canopy_closure, SR_leaf_litter_depth,
             SR_understory_height, SR_No.WaterBodies,
             ncol = 2)

## exploratory modelling ----

# linear model with all terms.  log transforming the Estimator may be justified
mod <- lm(log(diversity_estimate+1) ~ compartment_age + 
            canopy_closure + leaf_litter_depth + understory_height + No.WaterBodies,
            data = SR_data)

autoplot(mod, smooth.colour = NA)

# signficnant (almost) age and leaf litter.
# log transformation removes age.
summary(mod)

mod_step <- stepAIC(mod)
summary(mod_step)

etasq(mod, anova = TRUE)
etasq(mod_step, anova = TRUE)

# plot model results ----
# Step 1
# newX are defined on a grid of age and leaf litter (significant terms)
# and tat the mean of canopy and understory (insignficant)
newX <- expand.grid(compartment_age = seq(from = 0, to = 50, by = 1),
                    canopy_closure = mean(SR_data$canopy_closure),
                    leaf_litter_depth = mean(SR_data$leaf_litter_depth),
                    understory_height = seq(from = 30, to = 100, by = 10))

# Step 2: get the predictions
newY <- predict(mod, newdata = newX, interval = 'confidence')

# Step 3: Housekeeping to create data frame to plot
plotThese<-data.frame(newX, newY) %>% rename("Predicted Diversity" = fit)
head(plotThese)
# Step 4: make the plot
ggplot(plotThese, aes(x = compartment_age, y = `Predicted Diversity`, 
                      colour = factor(understory_height)))+
  geom_line() +
  theme_bw()
