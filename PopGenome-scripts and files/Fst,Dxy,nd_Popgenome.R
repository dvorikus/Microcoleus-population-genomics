#Calculate Fst, Dxy, and nucleotide diversity for dataset III with PopGenome

library(ape)
library(tidyverse)
library(PopGenome)
library(dplyr)
library(ggplot2)

#Load the vcf file with all strains (dataset III; without the outgroup)
microcoleus_vcf <- readVCF("./vcf/all.biallelic.version2.02.vcf.gz", numcols = 50000, tid = "CP003614.1", frompos = 1, topos = 7479014, include.unknown = TRUE)
populations_microcoleus_vcf <- read.csv("E:/Microcoleus PopGenome/populations_microcoleus_vcf.txt", sep="")
populations_vcf <- split(populations_microcoleus_vcf$ind, populations_microcoleus_vcf$pop)
microcoleus_vcf <- set.populations(microcoleus_vcf, populations_vcf, diploid = F)

#Transform our dataset to perform genome scans in sliding windows of 50kb with a step size of 12.5kb
microcoleus_sw_50 <- sliding.window.transform(microcoleus_vcf, width = 50000, jump = 12500, type = 2)

microcoleus_nuc_div_50 <- microcoleus_sw_50@nuc.diversity.within
microcoleus_sw_50_fst <- F_ST.stats(microcoleus_sw_50, mode = "nucleotide")
nd_50 <- microcoleus_sw_50_fst@nuc.diversity.within/50000
fst_value_50 <- t(microcoleus_sw_50_fst@nuc.F_ST.pairwise)
dxy_50 <- t(microcoleus_sw_50_fst@nuc.diversity.between/50000)
position <- seq(from = 1, to = 7479014, by = 12500) + 50000
position <- seq(from = 1, to = 7479014, by = 12500)
window_stop <- position + 50000
sum(window_stop > 7479014)
position <- position[which(window_stop < 7479014)]
window_stop <- window_stop[which(window_stop < 7479014)]
windows <- data.frame(start=position, stop=window_stop, mid= position + (window_stop - position)/2)
x <- colnames(fst_value_50)
x <- sub("/", "_", x)
colnames(fst_value_50) <- paste0(x, "_fst")
colnames(dxy_50) <- paste0(x, "_dxy")
microcoleus_data_50 <- as_tibble(data.frame(windows, nd_50, fst_value_50, dxy_50))
write.csv(microcoleus_data_50, "./microcoleus data 50.csv")

#Prepare separate files in excel to plot pairwise Fst and Dxy in a loop

fst_plots <- read_excel("fst_plots.xlsx")
colNames <- names(fst_plots)[2:79]
for (i in colNames) {
  plt <- ggplot(fst_plots, aes_string(x=mid/10^6, y=i)) +
    geom_line(size=0.7,colour="royalblue4") + geom_point(size=1, col="royalblue4") + ylim(0, 1) + theme_minimal() + xlab("Position (Mb)")
  print(plt)
  ggsave(plt, file=paste0("plot_",i,".pdf"))
  Sys.sleep(2)
}

dxy_plots <- read_excel("dxy_plots.xlsx")
colNames_dxy <- names(dxy_plots)[2:79]
for (i in colNames_dxy) {
  plt_dxy <- ggplot(dxy_plots, aes_string(x=mid/10^6, y=i)) +
    geom_line(size=0.7,colour="royalblue4") + geom_point(size=1, col="royalblue4") + ylim(0, 0.025) + theme_minimal() + xlab("Position (Mb)")
  print(plt_dxy)
  ggsave(plt_dxy, file=paste0("plot_",i,".pdf"))
  Sys.sleep(2)

  #Plot all Fst scatterplots in a single file; Remove M13 from the rest of analysis and for plotting, as it contains only one strain
  my_plots_fst <- list()
  fst_plots_without13 <- fst_plots %>% select(-contains("pop13"))
  colNames_without13 <- colnames(fst_plots_without13)[2:67]
  #then loop so you can get all plots
for (i in colNames_without13) {
           my_plots_fst[[i]] <- ggplot(fst_plots_without13, aes_string(x=mid/10^6, y=i)) +
                   geom_line(size=0.2,colour="royalblue4") + geom_point(size=0.2, col="royalblue4") + ylim(0, 1) + theme_minimal() + xlab("") + ylab("") + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
}
  
grid.arrange(grobs=my_plots_fst, ncol=12)
