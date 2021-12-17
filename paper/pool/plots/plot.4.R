
library(ggplot2)

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

data <- read.table("results.4.txt", head=T)

data$D <- as.factor(data$D)
data$S <- as.factor(data$S)

data$estimation=factor(data$estimation, levels=c("count","fit_count","fit_afl"))

data <- data[which((data$S!=10) & (data$D!=0.1) & (data$estimation %in% c("count","fit_count","fit_afl")) ),]

# a <- ggplot(data=data, aes(x=D, y=value, color=estimation)) + geom_point(position=position_dodge(0.5), size=4) + facet_grid(S ~ .)

a <- ggplot(data=data, aes(x=D, y=rmse, color=estimation)) + geom_point(position=position_dodge(0.5), size=4) + facet_grid(S ~ reference) + scale_colour_viridis_d()


a

ggsave("plot.4.png")


# ALEX, ignore what follows here....


names <- c()
for (i in unique(data$D)) {
        for (j in c("pop","sample",as.character(unique(data$est)))) {
                names<- rbind(names, cbind(i,j))
        }
}


sfs <- as.matrix(read.table("results.4.20.sfs.txt", head=F))

par(mfrow=c(2,2))
depths <- as.character(unique(names[,1]))
for (i in 1:length(depths)) {
	ind <- which(names[,1] == depths[i])
	ind <- ind[c(1:4,7)]
	barplot(sfs[ind,2:(ncol(sfs)-1)], beside=T, legend=names[ind,2], main=paste(depths[i],"X",sep="",collapse=""))
}

sfs <- as.matrix(read.table("results.4.50.sfs.txt", head=F))

par(mfrow=c(2,2))
depths <- as.character(unique(names[,1]))
for (i in 1:length(depths)) {
        ind <- which(names[,1] == depths[i])
        ind <- ind[c(1:4,7)]
        barplot(sfs[ind,2:(ncol(sfs)-80)], beside=T, legend=names[ind,2], main=paste(depths[i],"X",sep="",collapse=""))
}






