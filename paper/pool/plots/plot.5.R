
library(ggplot2)

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

data <- read.table("results.4.K.txt", head=T)

colnames(data)=c("S", "D", "estimation", "K" )

data <- data[which(data$estimation %in% c("count","fit_count","fit_afl") ),]

data$estimation <- factor(data$estimation, levels=c("fit_count","fit_afl"))

data$D <- factor(data$D, levels=c("1","2","5"))


a <- ggplot(data=data, aes(x=D, y=K, color=estimation)) + geom_point(position=position_dodge(0.5), size=4) + facet_grid(S ~ .) + scale_colour_viridis_d()



a


ggsave("plot.4.K.png")



