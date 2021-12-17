
cat("S", "D", "F", "estimation", "value", "\n")

for (S in c(20,50)) {

for (F in c(0.1, 0.2, 0.3, 0.4, 0.5)) {

for (D in c(0.5, 1, 2, 5)) {

	res <- read.table(paste("Results/",S,".",F,".",D,".out.gz",sep="",collapse=""))

	truth <- read.table(paste("Data/",S,".",F,".",D,".txt",sep="",collapse=""))

	# take only values with data
	truth <- truth[res[,2],]

	# pop maf
	pop <- truth[,5]
	pop[which(pop>0.5)] <- 1 - pop[which(pop>0.5)]

	# maf
	true <- truth[,(S+6)]
	true[which(true>0.5)] <- 1- true[which(true>0.5)]

	# MLE
	mle <- res[,11]

	for (i in 1:length(mle)) {
		cat(S,D,F,"SAF",true[i],"\n") # sample allele frequency, true
		cat(S,D,F,"MLE",mle[i],"\n")
	}

}
}
}




