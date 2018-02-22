## Nada Chapters 2-3 iNext Modelling ----

## libraries ----
library(iNEXT)
library(tidyverse)
library(gridExtra)
library(ggfortify)

## Step 0: data import ----

# forest compartment
compart <- read.csv('SpeciesAdultData/All_DTC2.csv')
habitat <- read.csv('SpeciesAdultData/DTC-AverageHabitatCompt-Age.csv')

str(compart)
str(habitat)

# Step 1: calculate column sums by forest category and compartment name to get iNEXT vals ----

# isolate the data on logged 1 time or unlogged
# get category and compartment name and species
# sort and summarise_all creates column sums
# age a new column that combines age and name (e.g. the three compartments within each age category)

comp_age<-compart %>% 
  filter(!str_detect(ForestType, "twice")) %>%
  select(ForestCategory, CompartmentName, Sp.01:Sp.27) %>%
  group_by(ForestCategory, CompartmentName) %>%
  summarise_all(.funs = function(x) (sum(x, na.rm = TRUE))) %>%
  mutate(Age_Name = paste(ForestCategory, CompartmentName, sep = ':'))

# Step 2: create Age based list for iNEXT ----

# for counting # 18 compartments
age_counter<-length(comp_age$CompartmentName)

# collection zone setup
compartment_input<- list()

# isolate the species count data for each age-compartment
for (i in 1:age_counter){
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
# some odd errors, as all sites have more than 1
compartment_mod <- iNEXT(compartment_use, datatype = 'abundance', nboot = 999)


# What we really want to do is to compute the estimate at a fixed min coverage

# at the minimum coverage
coverage_min_diversity <- estimateD(compartment_use, base = 'coverage') %>% 
  rename(Site = site) 

# at 80% coverage (extrapolated)
#coverage_min_diversity <- estimateD(compartment_use, base = 'coverage', level = 0.8) %>% 
#  rename(Site = site)

# visualise the results
# feel free to create more than one version of this!
par(mfrow = c(1,2))
plot(compartment_mod, type = 1)
plot(compartment_mod, type = 3)

# Step 4: Collect Asymptotic estimates for downstream analyses and add habitat data ----
diversity_data_obs <- compartment_mod$AsyEst %>% arrange(as.character(Site))
diversity_data_cov <- coverage_min_diversity %>% arrange(as.character(Site))
  
unique(diversity_data_cov$Site)
habitat2 <- arrange(habitat, Compartment) %>% rename(Site = Compartment)
habitat_use <- filter(habitat2, Site %in% diversity_data_cov$Site)
names(habitat_use)

# get everything together.  We need the habitat data represented three times each to match the 
# diversity data, which has the sites 3 times each for SR, Shannon and Simpson

df <- data.frame(diversity_data_cov, 
                 compartment_age = rep(habitat_use$ForestAge, each= 3),
                 canopy_closure = rep(habitat_use$CanopyClosure, each = 3),
                 leaf_litter_depth = rep(habitat_use$LeafLitterDepth, each = 3),
                 understory_height = rep(habitat_use$HerbPlantHeight, each = 3)) %>% 
  rename(coverage = SC) %>% rename(diversity_estimate = qD)

head(df)

## Exploratory Plotting of diversity and habitat data ----

SR_data <- filter(df, order == 0) %>% filter(compartment_age!='Unlogged') %>% 
  mutate(compartment_age = as.numeric(as.character(SR_data$compartment_age)))

SR_logged <- filter(df, order == 0) %>% filter(compartment_age =='Unlogged')

head(SR_data)
tail(SR_data)

# fascinating groups....
ggplot(SR_data, aes(x = compartment_age, y = diversity_estimate))+
  geom_point(size = 5)+
  geom_point(data = SR_logged, aes(x = 50, y = diversity_estimate), colour = 'red', size = 5)+
  geom_text(aes(label=Site),hjust=0.5, vjust=-2)+
  theme_bw()


## exploratory modelling ----

# linear model with all terms.  log transforming the Estimator may be justified
mod <- lm(log(diversity_estimate) ~ compartment_age + 
            canopy_closure + leaf_litter_depth + understory_height,
            data = SR_data)

autoplot(mod, smooth.colour = NA)

# signficnant (almost) age and leaf litter.
# log transformation removes age.
summary(mod)

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

