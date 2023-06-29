###### bactaxR code, ANI analyses

# instalation
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("ggtree", force = T)
# library(devtools)
# install_github("lmc297/bactaxR")

library(bactaxR)

ani <- read.ANI(file = "batch1to7.fastaANI.res")
summary(ani)

h <- ANI.histogram(bactaxRObject = ani, bindwidth = 0.1)
h

dend <- ANI.dendrogram(bactaxRObject = ani, ANI_threshold = 95, xline = c(4,5,6,7.5), xlinecol = c("#ffc425", "#f37735", "deeppink4", "black"), label_size = 0.5)

metadata <- dend$cluster_assignments$Cluster
names(metadata) <- dend$cluster_assignments$Genome

# color palette can be changed see: https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
ANI.graph(bactaxRObject = ani, ANI_threshold = 95,
          metadata = metadata,
          legend_pos_x = -1.5, show_legend = T, graphout_niter = 1000000, 
          legend_ncol = 1, color_palette =  rainbow(n = length(unique(metadata))), edge_color = "black")

write.table(metadata, file = "clusters.txt")