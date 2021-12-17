
# plot analysis 3

library(ggplot2)

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

data <- read.table("confused.txt", head=T)

data$F <- as.factor(data$F)
data$D <- as.factor(data$D)



a <- ggplot(data=data, aes(x=F, y=F1, color=D)) + geom_point(position=position_dodge(0.5), size=4) + facet_grid(. ~ S) + scale_colour_viridis_d()

a

ggsave("plot.3.png")







