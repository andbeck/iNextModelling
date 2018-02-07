## Nada Chapters 2-3 iNext Modelling ----

## libraries ----
library(iNEXT)
library(tidyverse)
library(gridExtra)

## Step 0: data import ----

# forest compartment
compart <- read.csv('SpeciesAdultData/All_DTC2.csv')

str(compart)

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

# for counting
age_counter<-length(unique(comp_age$Age_Name))

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
names(compartment_input)<-unique(comp_age$Age_Name)

# Get Rid of entries with SR < 1 (iNEXT does not work without >2 species)
compartment_input <- compartment_input[lapply(compartment_input, length)>1]

#  Step 3: iNEXT modelling happens here -----

# run iNEXT (it works on a list)
# some odd errors, as all sites have more than 1
compartment_mod <- iNEXT(compartment_input, datatype = 'abundance', nboot = 999)

# visualise the results
# feel free to create more than one version of this!
plot(compartment_mod)

# Step 4: Collect Asymptotic estimates for downstream analyses ----
# this is the data to which we will add habitat information.

# Collect Species Richness estimates (Shannon and Simpson also available)
SR_estimates <- compartment_mod$AsyEst %>%
  filter(Diversity == "Species richness")

# We used the combined coding of age and compartment
# Now we have to unsplit that to get an age only column back
# I call that new column Age_class
# str_split separates the text by the ':' in the Site variable produced by iNEXT
SR_estimates <- mutate(SR_estimates, 
                       Age_class = str_split(SR_estimates$Site, ":", simplify = TRUE)[,1])

# This creates mean and se estimates from the above replicated data on species richness
SR_sumDat <- SR_estimates %>%
  group_by(Age_class) %>%
  summarise(
    meanSR = mean(Estimator),
    seSR = sd(Estimator)/sqrt(n()))

# Step 5: Visualise the pattern of SR estimates with Age of forest since logging ----
p1 <- ggplot(SR_estimates, aes(x = Age_class, y = Estimator))+
  geom_point(size = 4, colour = 'cornflowerblue')+
  theme_bw()

p2 <- ggplot(SR_sumDat, aes(x = Age_class, y = meanSR, ymin = meanSR - seSR, ymax = meanSR + seSR))+
  geom_point(size = 4, colour = 'cornflowerblue')+
  geom_errorbar(width = 0.1)+
  theme_bw()

grid.arrange(p1, p2)

# Step 6: Add habitat data to the data frame (NADA TO DO) and then analyse using lm() and anova() ----

modSR <- lm(Estimator ~ Age_class, data = SR_estimates)
anova(modSR) # currently, the estimates do not vary by age class.