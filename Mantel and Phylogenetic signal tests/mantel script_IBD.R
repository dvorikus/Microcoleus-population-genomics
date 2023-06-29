#Script for testing the correlations between geographic and genetic distances of Microcoleus strains (Mantel)

library(geosphere)
library(vegan)
library(ggplot2)

#genetic distance matrix is the matrix of ANI values turned into binary matrix (0-1).

genetic.distance.matrix <- read.delim("genetic distance matrix.txt", header=FALSE)
coords <- read.csv("coords.csv")

geographic_distance <- as.dist(distm(coords[,2:3], fun = distGeo))
genetic_distance <- as.dist(1 - genetic.distance.matrix)

mantel(genetic_distance, geographic_distance, permutations = 9999, na.rm = TRUE)

#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#mantel(xdis = genetic_distance, ydis = geographic_distance, permutations = 9999,      na.rm = TRUE) 

#Mantel statistic r: 0.4412 
#      Significance: 1e-04

mantelplot <- ggplot(mat, aes(x = gg/1000, y = nn)) + geom_point(size=2, alpha=0.5, colour = "lightgrey") + labs(x = "Geographic distance (km)", y= "Genetic distance") + geom_smooth(method = "lm", colour="#009999", alpha=3) + theme(axis.text.x = element_text(face = "bold", colour="black", size = 12), axis.text.y = element_text(face = "bold", colour = "black", size = 12),  axis.title = element_text(face = "bold", size = 12, colour = "black"), panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"))
mantelplot + theme_minimal()



