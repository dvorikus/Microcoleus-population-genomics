#Script for testing the correlations between geographic and genetic distances of Microcoleus strains (Mantel)

library(geosphere)
library(vegan)
library(ggplot2)
library(jvamisc)

#genetic distance matrix is the matrix of ANI values turned into binary matrix (0-1).
ANI <- read_excel("ANI.xlsx", col_names = FALSE)
min(ANI, na.rm = T) #minimal ANI is 82.42258
ANI <- ANI[,-1]
#Transform lower triangle matrix into a symmetric one
ANI$...2 <- as.numeric(ANI$...2)
ANI$...292 <- as.numeric(ANI$...292)
ANI_n <- upper2full(t(ANI), diagval = 100)
ANI_b <- ANI_n / 100
genetic.distance.matrix <- as.dist(1 - ANI_b)

coords <- read.csv("coords.csv")
geographic_distance <- as.dist(distm(coords[,2:3], fun = distGeo))

mantel(genetic.distance.matrix, geographic_distance, permutations = 9999, na.rm = TRUE)

#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#mantel(xdis = genetic.distance.matrix, ydis = geographic_distance, permutations = 9999,      na.rm = TRUE) 

#Mantel statistic r: 0.4383 
#      Significance: 1e-04

nn <- as.vector(genetic.distance.matrix)
gg <- as.vector(geographic_distance)
mat <- data.frame(gg, nn)
mantelplot <- ggplot(mat, aes(x = gg/1000, y = nn)) + geom_point(size=2, alpha=0.5, colour = "lightgrey") + labs(x = "Geographic distance (km)", y= "Genetic distance") + geom_smooth(method = "lm", colour="#009999", alpha=3) + theme(axis.text.x = element_text(face = "bold", colour="black", size = 12), axis.text.y = element_text(face = "bold", colour = "black", size = 12),  axis.title = element_text(face = "bold", size = 12, colour = "black"), panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"))
mantelplot + theme_minimal()



