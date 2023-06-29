#Script for the presence of phylogenetic signal in climatic niche space (or habitat) of Microcoleus strains
#Specifically, we estimate here Pagel's lambda and Bloomberg's K.

library(geiger)
library(phytools)
library(ape)
library(picante)

#load tree and correct it for 0 branch lengths
tree <- read.newick("SpeciesTreeAlignment.fa.contree")
tree$edge.length <- tree$edge.length + 10^-7
t1 <- root(tree, outgroup = "M2_D5", resolve.root = TRUE)
is.rooted(t1)
#TRUE
is.binary(t1)
#TRUE

#rescaled tree for estimating Pagel's lambda (different from no phylogenetic signal)
t1.lambda.0 <- rescale(t1, "lambda", 0)
#rescaled tree for estimating Pagel's lambda (different from 1; or Brownian Motion model)
t1.lambda.1 <- rescale(t1, "lambda", 1)

#also we can use here brownie.lite as an alternative

#Load dataset with all variables and match the names of the tips to the dataframe
All.variables.mantel <- read.csv("All variables mantel.csv")
All_variables_phylo <- as.data.frame(All.variables.mantel)
label <- tree$tip.label
All_variables_phylo <- All_variables_phylo[match(label, All_variables_phylo$Strain),]

#correct the dataframe; turn NAs into 0 as the phylogenetic signal can not be calculated for missing values
All_variables_phylo[is.na(All_variables_phylo)] <- 0

#Here is the example of Bloomberg's K and Pagel's lambda for one variable (global fertilizer)
#Loop the script to get the values for the rest of variables

#Calculate Bloomberg's K

nfertilizer <- as.vector(t(All_variables_phylo$Fertilizer.global))
names(nfertilizer) <- tree$tip.label
K_nfert <- phylosignal(nfertilizer, t1)

#Calculate Pagel's lambda

#Significantly different from 1
nfert_t1 <- fitContinuous(t1, nfertilizer, model = "lambda")
nfert_t_0 <- fitContinuous(t1.lambda.0, nfertilizer)
nfert_t1_1 <- fitContinuous(t1.lambda.1, nfertilizer)
nfert_LR_dif1 <- 2*(nfert_t1$opt$lnL - nfert_t1_1$opt$lnL)
pchisq(nfert_LR_dif1, df=1, lower.tail = FALSE)

#Significantly different from 0

nfert_LR_dif0 <- 2*(nfert_t1$opt$lnL - nfert_t_0$opt$lnL)
pchisq(nfert_LR_dif0, df=1, lower.tail = FALSE)
