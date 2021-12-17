
for (S in c(10,20,50,100,1000)) {
	for (D in c(0.1,0.5,1,2,5,10)) {
		fin <- paste("Data/",S,".",D,".mpileup.gz", sep="", collapse="")
		depth <- as.numeric(read.table(fin)[,4])
		cat(S,D,length(depth),(S*2)*D,":",mean(depth),sd(depth),"\n")
	}
}




