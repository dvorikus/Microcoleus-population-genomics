#Calculate Fst, Dxy, and nucleotide diversity for dataset III with PopGenome

library(ape)
library(tidyverse)
library(PopGenome)
library(dplyr)
library(ggplot2)

#Load the vcf file with all strains (dataset III; without the outgroup)
microcoleus_vcf <- readVCF("./vcf/all.biallelic.version2.02.vcf.gz", numcols = 10000, tid = "CP003614.1", frompos = 1, topos = 7479014, include.unknown = TRUE)
populations_microcoleus_vcf <- read.csv("populations_microcoleus_vcf.txt", sep="")
populations_vcf <- split(populations_microcoleus_vcf$ind, populations_microcoleus_vcf$pop)
microcoleus_vcf <- set.populations(microcoleus_vcf, populations_vcf, diploid = F)

#Transform our dataset to perform genome scans in sliding windows of 10kb with a step size of 2.5kb
microcoleus_sw_10 <- sliding.window.transform(microcoleus_vcf, width = 10000, jump = 2500, type = 2)

#Calculate neutrality statistics
microcoleus_sw_10 <- neutrality.stats(microcoleus_sw_10)
tajima_10 <- microcoleus_sw_10@Tajima.D
write.csv(tajima_10, "./tajima 10.csv")
FusF_10 <- microcoleus_sw_10@Fu.Li.F
write.csv(FusF_10, "./fusF 10.csv")

#plot neutrality statistics in boxplots
position_10 <- seq(from = 1, to = 7479014, by = 2500)
window_stop_10 <- position_10 + 10000
sum(window_stop_10 > 7479014)
position_10 <- position_10[which(window_stop_10 < 7479014)]
window_stop_10 <- window_stop_10[which(window_stop_10 < 7479014)]
windows_10 <- data.frame(start=position_10, stop=window_stop_10, mid= position_10 + (window_stop_10 - position_10)/2)
tajima_10_for_plot <- as_tibble(data.frame(windows_10, tajima_10))
tajima_10_boxplot <- tajima_10_for_plot[, -c(1:3)]
tajima_10_boxplot_g <- tajima_10_boxplot %>% select(contains("pop")) %>% gather(key="species", value = "TajimasD")
tajima_10_boxplot_g$species <- as.character(tajima_10_boxplot_g$species)
tajima_10_boxplot_g$species <- factor(tajima_10_boxplot_g$species, levels = unique(tajima_10_boxplot_g$species))
ggplot(tajima_10_boxplot_g, aes(species, TajimasD)) + geom_boxplot(outlier.size = 1, alpha = 1) + theme_classic()

FusF_10_boxplot <- as.data.frame(FusF_10)
colnames(FusF_10_boxplot) <- c("pop1", "pop2", "pop3", "pop4", "pop5", "pop6", "pop7", "pop8", "pop9", "pop10", "pop11", "pop12")
fusF_boxplot_g <- FusF_10_boxplot %>% select(contains("pop")) %>% gather(key="species", value = "FusF")
fusF_boxplot_g$species <- as.character(fusF_boxplot_g$species)
fusF_boxplot_g$species <- factor(fusF_boxplot_g$species, levels = unique(fusF_boxplot_g$species))
ggplot(fusF_boxplot_g, aes(species, FusF)) + geom_boxplot(outlier.size = 1, alpha = 1) + theme_classic() + ylim(-3, 3)

#T test for each lineage to check if the values are significantly different from 0
t.test(tajima_10_for_plot$pop.1)
t.test(FusF_10_boxplot$pop.1)
...
