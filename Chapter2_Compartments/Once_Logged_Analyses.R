# This is the 1x analyses using the compiled data ===========================
# the data contain habitat variables (means and CV's), raw estiamtes of SR etc
# and iNEXT estimates from the compartments analysis script.

# 5 July 2018

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

# these data contain the habitat and iNEXT data
master_data <- read_csv("./Chapter2_Compartments/DataSources/All compartments_master_05.07.18.csv")
dim(master_data)
master_data$LoggingRotation
# get the  once logged and unlogged data
once <- filter(master_data, LoggingRotation == "once" | LoggingRotation == "unlogged")
glimpse(once)
dim(once)

## Exploratory Plotting of diversity and habitat data ----

# pairs correlations
ggpairs(select(once,Age_Latest, contains("_Mean")))

# Diversity metrics vs age since logging -----------
g1 <- ggplot(once, aes(x = Age_Latest, y = order_0_at_90))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=CompartmentName), point.padding = 0.5) +
  labs(x = "Compartment Age (Years Since Last Logging)",
       y = "Species Richness Estimate") +
  theme_bw()
g2 <- ggplot(once, aes(x = Age_Latest, y = order_1_at_90))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=CompartmentName), point.padding = 0.5) +
  labs(x = "Compartment Age (Years Since Last Logging)",
       y = "Shannon Richness Estimate") +
  theme_bw()
g3 <- ggplot(once, aes(x = Age_Latest, y = order_1_at_90))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=CompartmentName), point.padding = 0.5) +
  labs(x = "Compartment Age (Years Since Last Logging)",
       y = "Simpson Richness Estimate") +
  theme_bw()

grid.arrange(g1, g2, g3)


# Habitat Variable plots ----
SR_age <- ggplot(once, aes(x = Age_Latest, y = order_0_at_90))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=CompartmentName), point.padding = 0.5) +
  labs(x = "Compartment Age (Years Since Last Logging)",
       y = "Species Richness Estimate") +
  theme_bw()

SR_canopy_closure <- ggplot(once, aes(x = Canopy_Mean, y = order_0_at_90))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=CompartmentName), point.padding = 0.5) +
  labs(x = "Canopy Closure",
       y = "Species Richness Estimate") +
  theme_bw()

SR_leaf_litter_depth <- ggplot(once, aes(x = Leaflitter_Mean, y = order_0_at_90))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=CompartmentName), point.padding = 0.5) +
  labs(x = "Leaf Litter Depth",
       y = "Species Richness Estimate") +
  theme_bw()

SR_understory_height <- ggplot(once, aes(x = Herb_Mean, y = order_0_at_90))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=CompartmentName), point.padding = 0.5) +
  labs(x = "Understory Height",
       y = "Species Richness Estimate") +
  theme_bw()

SR_No_WaterBodies <- ggplot(once, aes(x = Water_Mean, y = order_0_at_90))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=CompartmentName), point.padding = 0.5) +
  labs(x = "No.WaterBodies",
       y = "Species Richness Estimate") +
  theme_bw()

SR_Primary <- ggplot(once, aes(x = NearPrimary, y = order_0_at_90))+
  geom_point(size = 5)+
  geom_text_repel(aes(label=CompartmentName), point.padding = 0.5) +
  labs(x = "near_Primary",
       y = "Species Richness Estimate") +
  theme_bw()

grid.arrange(SR_age, SR_canopy_closure, 
             SR_leaf_litter_depth, SR_understory_height, 
             SR_No_WaterBodies, SR_Primary,
             ncol = 3)

## linear modelling with iNEXT Asymptotes ----

## set the unlogged age: 100 or 200 ----
unlogged_age <- 200
once$Age_Latest <- iflese(once$Age_Latest == 200, unlogged_age, once$Age_Latest)

# Species Richness Analysis ---
# linear model with all terms; polynomial 2 and 1 
mod_full_poly <- lm(order_0_at_90 ~ poly(Age_Latest,2) + 
                      poly(Canopy_Mean,2) + poly(Leaflitter_Mean,2) + poly(Herb_Mean,2) + 
                      Water_Mean + NearPrimary, 
                    data = once)

mod_full_lin <- lm(order_0_at_90 ~ Age_Latest + 
                     Canopy_Mean + Leaflitter_Mean + Herb_Mean+ 
                     Water_Mean + NearPrimary, 
                   data = once)

# step 1 - no difference between models via SS/F
anova(mod_full_poly, mod_full_lin)

# step 2 - but the residuals of mod_full_lin are terrible
par(mfrow = c(2,2))
plot(mod_full_lin)

# step 3 - apply log() transformation to linear model
mod_full_lin_log <- lm(log(order_0_at_90) ~ Age_Latest + 
                         Canopy_Mean + Leaflitter_Mean + Herb_Mean+ 
                         Water_Mean + NearPrimary, 
                       data = once)

# step 4 diagnostics OK now
par(mfrow = c(2,2))
plot(mod_full_lin_log, add.smooth = FALSE)

# step 5: Anova - nothing is signficant
Anova(mod_full_lin_log)

# Step 6: Medium effect sizes....
eta <- etasq(mod_full_lin_log, anova = TRUE)
eta


#### START FIXING HERE! (5th July)

# Simpsons Diversity
# linear model with all terms; polynomial 2 and 1 
mod_full_poly_simp <- lm(diversity_estimate ~ poly(compartment_age,2) + 
                           poly(canopy_closure,2) + poly(leaf_litter_depth,2) + poly(understory_height,2) + 
                           No_WaterBodies + near_Primary, 
                         data = Simp_data)

mod_full_lin_simp <- lm(diversity_estimate ~ compartment_age + 
                          canopy_closure + leaf_litter_depth + understory_height + 
                          No_WaterBodies + near_Primary, data = Simp_data)

# step 1 - no difference between models via SS/F
anova(mod_full_poly_simp, mod_full_lin_simp)

# step 2: check diagnostics of simple model
par(mfrow = c(2,2))
plot(mod_full_lin_simp, add.smooth = FALSE)

# Step 3: Inference - the residials are fine without log for Simpsons
# and we find a singificant effect of compartment age!
# simpsons diversity increases with time since logging.
Anova(mod_full_lin_simp)
summary(mod_full_lin_simp)

# Step 4: effect sizes
# and time since logging has a large effect size too!
# and proximity to primary has a medium effect size.
eta_simp <- etasq(mod_full_lin_simp, anova = TRUE)
eta_simp


# Shannons Diversity
# linear model with all terms; polynomial 2 and 1 
mod_full_poly_shan <- lm(diversity_estimate ~ poly(compartment_age,2) + 
                           poly(canopy_closure,2) + poly(leaf_litter_depth,2) + poly(understory_height,2) + 
                           No_WaterBodies + near_Primary, 
                         data = Shan_data)

mod_full_lin_shan <- lm(diversity_estimate ~ compartment_age + 
                          canopy_closure + leaf_litter_depth + understory_height + 
                          No_WaterBodies + near_Primary, data = Shan_data)

# step 1 - no difference between models via SS/F
anova(mod_full_poly_shan, mod_full_lin_shan)

# step 2: check diagnostics of simple model
par(mfrow = c(2,2))
plot(mod_full_lin_shan, add.smooth = FALSE)

# Step 3: Inference - the residials are fine without log for Simpsons
# and we find a singificant effect of compartment age!
# simpsons diversity increases with time since logging.
Anova(mod_full_lin_shan)
summary(mod_full_lin_shan)

# Step 4: effect sizes
# and time since logging has a large effect size too!
# and proximity to primary has a medium effect size.
eta_shan <- etasq(mod_full_lin_shan, anova = TRUE)
eta_shan

effect_size_data <- data.frame(Term = rownames(eta), 
                               Species_Richness = eta$`Partial eta^2`,
                               Simpson = eta_simp$`Partial eta^2`,
                               Shannon = eta_shan$`Partial eta^2`)

effect_size_plot <- na.omit(
  gather(effect_size_data, key = q, value = Partial_eta, -Term))

ggplot(effect_size_plot, aes(x = Term, y = Partial_eta, colour = q))+
  geom_lollipop(size = 1)+
  ylab(expression(paste("Partial ", eta^2)))+
  facet_wrap(~q) +
  geom_hline(yintercept = c(0.01,0.06,0.28), linetype = 'dashed')+
  theme_bw(base_size = 15)+
  guides(colour = FALSE)+
  coord_flip()

out_anovas <- na.omit(data.frame(
  q = rep(c("Species Richness", "Shannon", "Simpson"), each = 7),
  rbind(
    tidy(eta), tidy(eta_shan), tidy(eta_simp))))

out_anovas <- out_anovas %>% mutate_if(is.numeric, funs(round(., 3))) %>% 
  rename(Term = term, `F` = statistic, SS = sumsq)

write_csv(out_anovas, path = "./Chapter2_Compartments/DataSources/out_anovas.csv")

# MVABUND approach ----

# does community composition vary by habitat.

# prepare the abundance data
abund_dat <- comp_age %>% 
  ungroup() %>% 
  select(starts_with("Sp.")) %>% 
  mvabund()

# prepare the habitat data
hab_dat <- SR_data %>% 
  select(-Site, -order, -coverage, - diversity_estimate) %>% 
  as.data.frame()

# view some of the abundance data
plot(abund_dat, hab_dat$near_Primary)
plot(abund_dat, cut(hab_dat$compartment_age, 3))

# build the multivariate model of site x species abundance ~ habitat
mv_model <- manyglm(abund_dat ~ compartment_age+
                      canopy_closure +
                      leaf_litter_depth +
                      understory_height +
                      No_WaterBodies +
                      near_Primary, 
                    data = hab_dat, family = "negbinomial")

# diagnostics (some hint of mean-var with poisson.... negBin much better)
par(mfrow = c(1,3))
plot(mv_model, which = 1:3)

# model inference
set.seed(260518)
anova(mv_model)

# betapart to dig deeper?
library(betapart)
library(vegan)

beta_dat <- comp_age %>% 
  ungroup() %>% 
  select(starts_with("Sp."))

# create beta dissimilarity based on abundance
betaCore <- beta.pair.abund(beta_dat)
betaCore <- bray.part(beta_dat)
betaCore

# the MVabund model highlighted leaf litter and No.Water Bodies
# lets create a 3 level factor for each of these, and assess beta
timeSince_3fac <- cut(hab_dat$compartment_age, breaks = c(0,18,50,200), 
                      labels = c("recent", "recovered", "primary-like")) 
canopy_3fac <- cut(hab_dat$canopy_closure, 3, 
                   labels = c("closed", "speckled", "open"))
litter_3fac <- cut(hab_dat$leaf_litter_depth, 3, 
                   labels = c("shallow", "med", "deep"))
understory_3fac <- cut(hab_dat$understory_height, breaks = c(0,65,80,100),
                       labels = c("low", "med", "tall"))
water_3fac <- cut(hab_dat$No_WaterBodies, 3, 
                  labels = c("few", "mid", "many"))

par(mfrow = c(2,3))
plot(betadisper(betaCore[[3]], timeSince_3fac), main = "Time Since Logging")
plot(betadisper(betaCore[[3]], canopy_3fac), main = "Canopy Closure")
plot(betadisper(betaCore[[3]], litter_3fac), main = "Litter Depth **")
plot(betadisper(betaCore[[3]], understory_3fac), main = "Understory Height")
plot(betadisper(betaCore[[3]], water_3fac), main = "Water Bodies **")
plot(betadisper(betaCore[[3]], hab_dat$near_Primary), main = "Near Primary")

# betapart approach does not show any signficant changes
anova(betadisper(betaCore[[3]], timeSince_3fac))
anova(betadisper(betaCore[[3]], canopy_3fac))
anova(betadisper(betaCore[[3]], litter_3fac))
anova(betadisper(betaCore[[3]], understory_3fac))
anova(betadisper(betaCore[[3]], water_3fac))
anova(betadisper(betaCore[[3]], hab_dat$near_Primary))


# # model averaging approach ----
# library(MuMIn)
# 
# fullMod_dredge <- update(mod_full_lin, na.action = na.fail)
# dd <- dredge(fullMod_dredge, subset = 
#                dc(compartment_age, I(compartment_age^2)) && 
#                dc(canopy_closure, I(canopy_closure^2)) &&
#                dc(leaf_litter_depth, I(leaf_litter_depth^2)) && 
#                dc(understory_height, I(understory_height^2)))
#   
# 
# # set for delta AICc < 2 : the AICc is compared to the best model
# # we want ones that are close to the best model
# # delta = 'difference in'
# summary(model.avg(dd, subset = delta < 3))
# 
# # best model
# summary(get.models(dd, 1)[[1]])


# # evenness ----
# ## Species richness (S) and Pielou's evenness (J):
# # S <- specnumber(BCI) ## rowSums(BCI > 0) does the same...
# # J <- H/log(S)
# 
# even <- function(x){
#   H <- diversity(x)
#   S <- sum(x>0)
#   J <- H/log(S)
#   return(J)
# }
# 
# # the pieces
# compartment_use$PB87
# diversity(compartment_use$PB87)
# sum(compartment_use$PB87>0)
# 
# # test on single compartment
# even(compartment_use$PB87)
# 
# # run on all compartments
# evens <- t(map_df(compartment_input, function(x) even(x))) %>% data.frame()
# 
# # craft data frame for plotting
# names(evens) <- "evenness"
# evens <- rownames_to_column(evens)
# names(evens)[1] <- "Compartment"
# evens <- arrange(evens, Compartment)
# habs <- arrange(habitat, Compartment)
# even_habs <- data.frame(evenness = evens$evenness, habs)
# 
# e1 <- ggplot(even_habs, aes(x = Age.of.forest, y = evenness, label = Compartment))+
#   geom_point(size = 3)+
#   geom_text_repel()
# 
# e2 <- ggplot(even_habs, aes(x = herb_mean_height, y = evenness, label = Compartment))+
#   geom_point(size = 3)+
#   geom_text_repel()
# 
# e3 <- ggplot(even_habs, aes(x = litter_mean_depth, y = evenness, label = Compartment))+
#   geom_point(size = 3)+
#   geom_text_repel()
# 
# e4 <- ggplot(even_habs, aes(x = canopy_mean_closure, y = evenness, label = Compartment))+
#   geom_point(size = 3)+
#   geom_text_repel()
# 
# e5 <- ggplot(even_habs, aes(x = water_mean_number, y = evenness, label = Compartment))+
#   geom_point(size = 3)+
#   geom_text_repel()
# 
# e6<- ggplot(even_habs, aes(x = AdjacentUnlogged, y = evenness, label = Compartment))+
#   geom_point(size = 3)+
#   geom_text_repel()
# 
# grid.arrange(e1,e2,e3,e4,e5,e6, ncol = 3)
