

library(ggplot2)


data <- read.table("results.4.K.txt", head=T)
colnames(data)=c("S", "D", "estimation", "K" )
data <- data[which(data$estimation %in% c("count","fit_count","fit_afl") ),]
data$estimation <- factor(data$estimation, levels=c("fit_count","fit_afl"))

#data$D <- factor(data$D, levels=c("1","2","5"))
labels_plot<-breaks_plot
labels_plot[!(breaks_plot%in%unique(data$D))]<-""

a <- ggplot(data=data, aes(x=D, y=K, color=estimation, group=estimation)) + 
  geom_line(size=0.5) + 
  #scale_shape_manual(values=c(19,21))+
  geom_point(size=2,fill="white",stroke=1) + 
  scale_colour_manual(values=c("#003F91","#58A4B0"))+
  scale_x_continuous(breaks=breaks_plot,labels=labels_plot,trans="log2")+
  labs(x="Depth",y="K")+
  facet_grid(. ~ S) +
  theme_minimal()+
  theme(panel.grid.major=element_line(color="gray95"),panel.grid.minor=element_line(color="gray95"))




