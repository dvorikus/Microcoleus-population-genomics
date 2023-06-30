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