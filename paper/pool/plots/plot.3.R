
# plot analysis 3

library(ggplot2)


data <- read.table("confused.txt", head=T)
data$D <- as.factor(data$D)

a <- ggplot(data=data, aes(x=F, y=F1, color=D,group=D)) + 
  geom_line(size=0.5) + 
  #scale_linetype_manual(values=c("dotdash", "dotted")) +
  geom_point(shape=19,size=2,fill="white",stroke=1) + 
  labs(x="Fixed population allele frequency (F)",y="F1")+
  facet_grid(. ~ S) + 
  scale_colour_viridis_d()+
  theme_minimal()+
  theme(panel.grid.major=element_line(color="gray95"),panel.grid.minor=element_line(color="gray95"))





