
# plot analysis 7

library(ggplot2)

data <- read.table("results.7.txt", head=T, stringsAsFact=F)
data[,2][which(data[,2]=="null")]<-"non_causal"
colnames(data)=c("D","SNP","LRT")
data$D=as.factor(data$D)
data$SNP=factor(data$SNP, levels=c("non_causal","causal"))
#data$LRT=as.factor(data$LRT)

a <- ggplot(data=data,aes(x=SNP, y=LRT, color=SNP,fill=SNP)) +
  geom_boxplot(position=position_dodge(0.85),width=1,alpha=0.5) + 
  labs(x="SNP",y="LRT")+
  scale_colour_manual(values=c("#4C4C4C","#DDD92A"))+
  scale_fill_manual(values=c("#4C4C4C","#DDD92A"))+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major=element_line(color="gray95"),panel.grid.minor=element_line(color="gray95"))+
  facet_grid(. ~ D)





