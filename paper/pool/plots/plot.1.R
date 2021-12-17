
# plot analysis 1

library(ggplot2)

# The palette with grey:
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

# To use for line and point colors, add scale_colour_manual(values=cbPalette)


data <- read.table("results.1.txt", head=T)

data2 <- data[which(data$S!=10 & data$D!=10 & data$D!=0.1),]

data2$D <- as.factor(data2$D)

a <- ggplot(data=data2, aes(x=D, y=value, shape=estimation, color=reference)) + geom_point(position=position_dodge(0.5), size=4) + facet_grid(error ~ S, scales="free") + scale_colour_viridis_d()

ggsave("plot.1.png")







