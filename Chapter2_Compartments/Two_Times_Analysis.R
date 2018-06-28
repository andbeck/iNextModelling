# 2x logging analysis

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

# data ----
two_times_habitat <- read_csv("./Chapter2_Compartments/DataSources/TwoTimes_iNEXT90_With Habitat Variables.csv")

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

grid.arrange(g1,g2,g3,g4)
