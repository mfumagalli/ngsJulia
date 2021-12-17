
# plot analysis 7

library(ggplot2)

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

data <- read.table("results.7.txt", head=T, stringsAsFact=F)

data[,2][which(data[,2]=="null")]<-"non_causal"

colnames(data)=c("D","SNP","LRT")

data$D=as.factor(data$D)
data$SNP=factor(data$SNP, levels=c("non_causal","causal"))
#data$LRT=as.factor(data$LRT)


a <- ggplot(data=data) + geom_boxplot(aes(x=SNP, y=LRT), notch=T) + facet_grid(. ~ D)

a

ggsave("plot.7.png")







