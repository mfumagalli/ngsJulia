
# plot analysis 1

library(ggplot2)


data <- read.table("results.1.txt", head=T)
data2 <- data[which(data$S!=10 & data$D!=10 & data$D!=0.1),]

breaks_plot<- seq(0,5.5,0.5) # 2^(seq(-1,2.5,0.5))
labels_plot<-breaks_plot
labels_plot[!(breaks_plot%in%unique(data2$D))]<-""

a<-ggplot(data=data2, aes(x=D, y=value, shape=reference,linetype=reference,group=paste(reference,estimation),color=estimation)) + 
  geom_line(size=0.5) + 
  #scale_linetype_manual(values=c("dotdash", "dotted")) +
  scale_shape_manual(values=c(19,21))+
  geom_point(size=2,fill="white",stroke=1) + 
  facet_grid( error ~  S + estimation, scales="free")+
  scale_colour_manual(values=c("#D81159","#218380","#FFBC42"))+
  scale_x_continuous(breaks=breaks_plot,labels=labels_plot,trans="log2")+
  theme_minimal()+
  labs(x="Depth",y=NULL)+
  theme(panel.grid.major=element_line(color="gray95"),panel.grid.minor=element_line(color="gray95"))








