#Get correlations between the number of CRISPRs per strain, the r/m values, and the number of recombination blocks per strain (from Gubbins output)

CRISPRs_per_strains <- read_excel("CRISPRs per strains.xlsx")
crispr <- CRISPRs_per_strains$`Number of CRISPRs`

gubbins_stats <- read_excel("Supplementary Table S10.xlsx")
rm <- gubbins_stats$`r/m`

cor.test(crispr, rm)

#Pearson's product-moment correlation

#data:  crispr and rm
#t = 2.2475, df = 199, p-value = 0.02571
#alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.01936362 0.28942923
#sample estimates:
#      cor 
# 0.1573367 

number_recomb_blocks <- gubbins_stats$Num_of_Recombination_Blocks
cor.test(crispr, number_recomb_blocks)

#Pearson's product-moment correlation

#data:  crispr and number_recomb_blocks
#t = -2.8818, df = 199, p-value = 0.004388
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.32942056 -0.06351537
#sample estimates:
#       cor 
#-0.2001507 