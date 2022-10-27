
# plot analysis 2

library(ggplot2)


data <- read.table("results.2.txt", head=T)

data$F <- as.factor(data$F)

a <- ggplot(data=data,aes(x=F, y=value, color=estimation,fill=estimation)) + 
  geom_violin(position=position_dodge(0.85),trim=FALSE,adjust=4,alpha=0.5,color=NA)+
  geom_boxplot(position=position_dodge(0.85),fill="white",width=0.3,outlier.shape=NA) + 
  facet_grid(S ~ D) +
  labs(x="Fixed population allele frequency (F)",y="Estimated minor allele frequency")+
  scale_colour_manual(values=c("#218380","#503F4E"))+
  scale_fill_manual(values=c("#218380","#503F4E"))+
  theme_minimal()+
  theme(panel.grid.major=element_line(color="gray95"),panel.grid.minor=element_line(color="gray95"))





