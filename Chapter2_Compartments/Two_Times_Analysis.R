# 2x logging analysis

library(heplots) # effect sizes +
library(mvabund) # 
library(vegan) #

library(gridExtra) # multi-panel ggplots
library(ggfortify) # diagnostics for models
library(ggrepel) # nice labelling of points
library(GGally)
library(broom)
library(tidyverse)

# data ----
two_times_habitat <- read_csv("./Chapter2_Compartments/DataSources/All compartments_master_05.07.18.csv")

glimpse(two_times_habitat)
two_times_habitat$Age_Latest

# initial graphs
g1 <- ggplot(two_times_habitat, aes(x = LoggingRotation, y = SR_per_compartment))+
  geom_boxplot()

g2 <- ggplot(two_times_habitat, aes(x = Age_Latest, y = SR_per_compartment, 
                              colour = LoggingRotation,
                              label = CompartmentName))+
  geom_point()+
  geom_label_repel()+
  facet_wrap(~LoggingRotation)+
  theme_bw(base_size = 15)+
  theme(legend.position = 'none')


g3 <- ggplot(two_times_habitat, aes(x = LoggingRotation, y = order_0_at_90))+
  geom_boxplot()

g4 <- ggplot(two_times_habitat, aes(x = Age_Latest, y = order_0_at_90, 
                              colour = LoggingRotation,
                              label = CompartmentName))+
  geom_point()+
  geom_label_repel()+
  facet_wrap(~LoggingRotation)+
  theme_bw(base_size = 15)+
  theme(legend.position = 'none')

g5 <- ggplot(two_times_habitat, aes(x = LoggingRotation, y = SR_Asymptote))+
  geom_boxplot()

g6 <- ggplot(two_times_habitat, aes(x = Age_Latest, y = SR_Asymptote, 
                                    colour = LoggingRotation,
                                    label = CompartmentName))+
  geom_point()+
  geom_label_repel()+
  facet_wrap(~LoggingRotation)+
  theme_bw(base_size = 15)+
  theme(legend.position = 'none')


grid.arrange(g1, g2, g3, g4, g5, g6, nrow = 3)

# Formal Analysis

# get the data for once/twice and age <=20
comp_1x_2x <- filter(two_times_habitat, LoggingRotation == 'once' | LoggingRotation == "twice") %>% 
  filter(Age_Latest<=20)

glimpse(comp_1x_2x)
comp_1x_2x$order_0_at_90

# plot it (choose metric!)
ggplot(comp_1x_2x, aes(x = Age_Latest, y = SR_per_compartment, colour = LoggingRotation))+
  geom_point(size = 5)+
  geom_smooth(method = lm)+
  theme_bw()

# model it/
model <- lm(SR_per_compartment ~ Age_Latest * LoggingRotation, data = comp_1x_2x)
Anova(model)
summary(model)

