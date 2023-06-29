# install
# devtools::install_github("gtonkinhill/fastbaps")
# The fast BAPS algorithm is based on applying the hierarchical Bayesian clustering (BHC) algorithm of [@Heller2005-kp] to the problem of clustering genetic sequences using the same likelihood as BAPS [@Cheng2013-mp]. The Bayesian hierarchical clustering can be initiated with sequences as individual clusters or by running a faster conventional hierarchical clustering initially followed by BHC of the resulting clusters.
 
# The algorithm has been written to take advantage of fast sparse matrix libraries and is able to handle 1000's of sequences and 100,000's of SNPs in under an hour on a laptop using a single core.

# Alternatively, we can condition on an initial phylogentic or hierarchical tree and provide the partition of the hierarchy that maximises the BAPS likelihood. This is useful if the user is mainly interested in partitioning an already calculated phylogeny. We have also noticed that partitioning a hierarchy built using ward.D2 distance gives very reasonable results, very quickly.

# Libraries

library(fastbaps)
library(ggtree)
library(phytools)
library(ggplot2)


## Loading data

# We first need to load a multiple sequence alignment into sparse format. We can choose between the original BAPS prior or a prior proportional to the mean frequency of each allele in the population.

sparse.data <- import_fasta_sparse_nt("SpeciesTreeAlignment.fa")


# Here we make use of the 'optimised symmetric' prior, which empirically chooses the variance of the Dirichlet prior on the component mixtures.


sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric")


## Running fastbaps

#It is a good idea to choose `k.init` to be significantly larger than the number of clusters you expect. By default it is set to the number of sequences / 4.


baps.hc <- fast_baps(sparse.data, k.init=80)

# This provides a Bayesian hierarchical clustering of the data. To obtain the partition of this hierarchy under Dirichlet Process Mixture model run

best.partition <- best_baps_partition(sparse.data, baps.hc)


#We can  plot the output of the  algorithm along with a pre-calculated tree using ggtree [@Yu2017-bf].


iqtree <- phytools::read.newick("SpeciesTreeAlignment.fa.treefile")
plot.df <- data.frame(id=colnames(sparse.data$snp.matrix),
                      fastbaps=best.partition,
                      stringsAsFactors = FALSE)

gg <- ggtree(iqtree)

f2 <- facet_plot(gg, panel="fastbaps", data=plot.df, geom=geom_tile, aes(x=fastbaps), color='blue')
f2


# We can compare this result to other priors, the un-optimised symmetric or BAPS prior similar to STRUCTURE and hierBAPS, an optimised BAPS prior or the population mean based prior of Heller et al.

sparse.data <- optimise_prior(sparse.data, type = "baps")

baps.hc <- fast_baps(sparse.data)
best.partition <- best_baps_partition(sparse.data, baps.hc)

plot.df <- data.frame(id=colnames(sparse.data$snp.matrix),
                      fastbaps=best.partition,
                      stringsAsFactors = FALSE)

gg <- ggtree(iqtree)
f2 <- facet_plot(gg, panel="fastbaps", data=plot.df, geom=geom_tile, aes(x=fastbaps), color='blue')
f2


# we can also use the same prior as used in the BHC algorithm of Heller et al. However this tends to overpartition population genetic data.


sparse.data <- optimise_prior(sparse.data, type = "hc")

baps.hc <- fast_baps(sparse.data)
best.partition <- best_baps_partition(sparse.data, baps.hc)

plot.df <- data.frame(id=colnames(sparse.data$snp.matrix),
                      fastbaps=best.partition,
                      stringsAsFactors = FALSE)

gg <- ggtree(iqtree)
f2 <- facet_plot(gg, panel="fastbaps", data=plot.df, geom=geom_tile, aes(x=fastbaps), color='blue')
f2


# we can also investigate multiple levels


sparse.data <- import_fasta_sparse_nt("SpeciesTreeAlignment.fa")
multi <- multi_res_baps(sparse.data)

plot.df <- data.frame(id=colnames(sparse.data$snp.matrix),
                      fastbaps=multi$`Level 1`,
                      fastbaps2=multi$`Level 2`,
                      stringsAsFactors = FALSE)

gg <- ggtree(iqtree)

f2 <- facet_plot(gg, panel="fastbaps level 1", data=plot.df, geom=geom_tile, aes(x=fastbaps), color='blue')
f2 <- facet_plot(f2, panel="fastbaps level 2", data=plot.df, geom=geom_tile, aes(x=fastbaps2), color='green')
f2


# We can also partition an initial hierarchy or phylogeny.


sparse.data <- import_fasta_sparse_nt("SpeciesTreeAlignment.fa", prior = "baps")

iqtree.rooted <- phytools::midpoint.root(iqtree)
best.partition <- best_baps_partition(sparse.data, iqtree.rooted)

plot.df <- data.frame(id=iqtree.rooted$tip.label,
                      fastbaps=best.partition,
                      stringsAsFactors = FALSE)

gg <- ggtree(iqtree.rooted)
f2 <- facet_plot(gg, panel="fastbaps", data=plot.df, geom=geom_tile, aes(x=fastbaps), color='blue')
f2

# We can compare this result to other priors, the un-optimised symmetric or BAPS prior similar to STRUCTURE and hierBAPS, an optimised BAPS prior or the population mean based prior of Heller et al.
# with mid point
# there two dataset for the best.partition - sparse.data, baps.hc (it is from the the optimaztion step above)
sparse.data <- optimise_prior(sparse.data, type = "baps")

baps.hc <- fast_baps(sparse.data)
best.partition <- best_baps_partition(sparse.data, baps.hc)

plot.df <- data.frame(id=colnames(sparse.data$snp.matrix),
                      fastbaps=best.partition,
                      stringsAsFactors = FALSE)

gg <- ggtree(iqtree.rooted)
f2 <- facet_plot(gg, panel="fastbaps", data=plot.df, geom=geom_tile, aes(x=fastbaps), color='blue')
f2



# finally we can also look at the stability of the inferred clusters using the Bootstrap

sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric")
boot.result <- boot_fast_baps(sparse.data)
dendro <- as.dendrogram(fast_baps(sparse.data))
gplots::heatmap.2(boot.result, dendro, dendro, tracecol=NA)

write.table(plot.df, file="fastbaps.populations.unoptimazed.priors.similar.to.Structure.heirBAPS.result.txt")
## References

sessionInfo()
