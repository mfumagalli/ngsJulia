
#args <- commandArgs(T)
#fin <- args[1]
#fout <- args[2]
#rm(args)

fout="test.likes.pdf"

pdf(file=fout)

for (i in c(2,5,10,20,50)) {

	fin=paste("Results/test.",i,".likes",sep="",collapse="")

	res <- read.table(fin, sep="\t", head=T, stringsAsFact=F)

	midPoints <- res[,1]+(res[,2]-res[,1]+1)/2

	true <- rep(c(1,2,3,4,5,5,4,3,2,1), each=10)

	plot(x=midPoints, y=true, type="s", ylab="Ploidy", xlab="Genome", col="red", cex=2, main=paste("Depth:",i), sub="MLE")

	points(x=midPoints, y=res$ploidy, ty="s", col="black", lty=2, cex=2)

	legend("topright", col=c("red","black"), lty=c(1,2), legend=c("True","Est."))

}

dev.off()

fout="test.pprobs.pdf"

pdf(file=fout)

for (i in c(2,5,10,20,50)) {

	fin=paste("Results/test.",i,".pprobs",sep="",collapse="")

	res <- read.table(fin, sep="\t", head=T, stringsAsFact=F)

	midPoints <- res[,1]+(res[,2]-res[,1]+1)/2

	true <- rep(c(1,2,3,4,5,5,4,3,2,1), each=10)

	plot(x=midPoints, y=true, type="s", ylab="Ploidy", xlab="Genome", col="red", cex=2, main=paste("Depth:",i), sub="Bayes")

	points(x=midPoints, y=res$ploidy, ty="s", col="black", lty=2, cex=2)

	legend("topright", col=c("red","black"), lty=c(1,2), legend=c("True","Est."))

}

dev.off()


