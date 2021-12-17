
# plot analysis 2

library(ggplot2)

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)


data <- read.table("results.2.txt", head=T)

data$F <- as.factor(data$F)

a <- ggplot(data=data) + geom_boxplot(aes(x=F, y=value, fill=estimation), notch=T) + facet_grid(D ~ S)

a

ggsave("plot.2.png")







