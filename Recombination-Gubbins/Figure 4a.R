#Making a bubble plot with ggplot2

library(ggplot2)
library(reshape2)

#load csv file with recombination fractions between populations
between_species <- read.csv("genome fractions - bubble plot")

#melt the data, so it will be a data frame for plotting
melted_between_species <- melt(between_species)

#change both variable (columns where the names of the populations are) into characters
#so they will appear in the order in the plot
melted_between_species$X <- as.character(melted_between_species$X)
melted_between_species$variable <- as.character(melted_between_species$variable)

melted_between_species$X <- factor(melted_between_species$X, levels = unique(melted_between_species$X))
melted_between_species$variable <- factor(melted_between_species$variable, levels = unique(melted_between_species$variable))

#then we can plot

sizeRange <- c(0, 10) #size of the bubble
ggplot(melted_between_species, aes(x = X, y = variable)) + geom_point(aes(size = value, colour = value))+ scale_colour_gradient(low = "green", high = "red") + theme_minimal() + scale_size(range = sizeRange)
