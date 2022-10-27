
library(ggplot2)

data <- read.table("results.4.txt", head=T)

#data$D <- as.factor(data$D)
data$S <- as.factor(data$S)

data$estimation=factor(data$estimation, levels=c("count","fit_count","fit_afl"))

data <- data[which((data$S!=10) & (data$D!=0.1) & (data$estimation %in% c("count","fit_count","fit_afl")) ),]

# a <- ggplot(data=data, aes(x=D, y=value, color=estimation)) + geom_point(position=position_dodge(0.5), size=4) + facet_grid(S ~ .)
breaks_plot<- seq(0,5.5,0.5) # 2^(seq(-1,2.5,0.5))
labels_plot<-breaks_plot
labels_plot[!(breaks_plot%in%unique(data$D))]<-""

a <- ggplot(data=data, aes(x=D, y=rmse, linetype=reference,shape=reference,color=estimation, group=paste(reference,estimation))) + 
  geom_line(size=0.5) + 
  scale_shape_manual(values=c(19,21))+
  geom_point(size=2,fill="white",stroke=1) + 
  scale_colour_manual(values=c("#F7B801","#003F91","#58A4B0"))+
  scale_x_continuous(breaks=breaks_plot,labels=labels_plot,trans="log2")+
  labs(x="Depth",y="RMSE")+
  facet_grid( reference  ~  S, scales="free")+
  theme_minimal()+
  theme(panel.grid.major=element_line(color="gray95"),panel.grid.minor=element_line(color="gray95"))






