#Script for plotting r/m and rho/theta values from the Gubbins output
#Plot rho/theta values
library(tidyverse)
library(dplyr)
library(ggplot2)
library(Hmisc)

rho_theta_per_pop <- read_excel("rho_theta_per_pop.xlsx")
ro_t_g <- rho_theta_per_pop %>% select(contains("pop")) %>% gather(key="species", value="pop")

#for ordering them in the correct order from pop1-pop12
#turn your species column into a character vector
ro_t_g$species <- as.character(ro_t_g$species)
#then turn it back into a factor with the levels in the correct order
ro_t_g$species <- factor(ro_t_g$species, levels = unique(ro_t_g$species))

#then plot
ggplot(ro_t_g, aes(species, pop)) + geom_boxplot(outlier.size = 1, alpha = 1) + theme_minimal() + ylab("Rho/theta") + xlab("")

#Now for the r/m
rm_per_pop <- read_excel("rm_per_pop.xlsx")
rm_t_g <- rm_per_pop %>% select(contains("pop")) %>% gather(key="species", value = "pop")
rm_t_g$species <- as.character(rm_t_g$species)
rm_t_g$species <- factor(rm_t_g$species, levels = unique(rm_t_g$species))
ggplot(rm_t_g, aes(species, pop)) + geom_boxplot(outlier.size = 1, alpha = 1) + theme_minimal() + ylab("r/m") + xlab("")

#In case you want to get violin plots do this; trim is trimming peaks
testing <- ggplot(ro_t_g, aes(species, pop)) + geom_violin(trim = FALSE) + theme_minimal() + ylab("Rho/theta") + xlab("")
testing + stat_summary(fun.y = mean, geom = "point", color = "red", size=1)
#that is rho theta. Now down is for rm the same
testing_rm <- ggplot(rm_t_g, aes(species, pop)) + geom_violin(trim = FALSE) + theme_minimal() + ylab("r/m") + xlab("")
testing_rm + stat_summary(fun.y = mean, geom = "point", color = "red", size=1)
