
cat("D", "test", "LRT", "\n")

for (D in c(0.5, 1, 2, 5)) {

	for (F in c("null","causal")) {

		safs <- list()
		i <- 0

		for (S in c("cases","controls")) {

			i <- i+1
			safs[[i]] <- read.table(paste("Results/",S,".",F,".",D,".saf.gz",sep="",collapse=""))

		}

		p0 <- apply(FUN=which.max, X=safs[[2]], MAR=1) # in practice should be -1, these are indexs

		for (j in 1:length(p0)) {

			H0 <- safs[[1]][j,p0[j]]
			HA <- max(safs[[1]][j,])
			LRT <- -2*(H0-HA)

			cat(D,F,LRT,"\n")

		}
	}
}



