
# calculate estimation error

cat("S", "D", "sites", "error", "reference", "estimation", "value", "\n")

for (S in c(10,20,50,100)) { # nr samples

for (D in c(0.1, 0.5, 1, 2, 5, 10)) { # depth

	# estimates
	res <- read.table(paste("Results/",S,".",D,".out.gz",sep="",collapse=""))

	# true values
	truth <- read.table(paste("Data/",S,".",D,".txt",sep="",collapse=""))
	# take only values with data, therefore with estimates
	truth <- truth[res[,2],]

	# nr valid sites
	nsites <- nrow(truth)

	# true population maf
	pop <- truth[,5]
	pop[which(pop>0.5)] <- 1 - pop[which(pop>0.5)]

	# true sample maf
	true <- truth[,(S+6)]
	true[which(true>0.5)] <- 1- true[which(true>0.5)]

	# rmse against pop and sample values
	rmse_maf <- c(mean(sqrt((true-res[,10])^2)), mean(sqrt((true-res[,11])^2)), mean(sqrt((true-res[,12])^2)))
	rmse_maf_pop <- c(mean(sqrt((pop-res[,10])^2)), mean(sqrt((pop-res[,11])^2)), mean(sqrt((pop-res[,12])^2)))

	# bias pred-true
	sb_maf <- c(mean(-(true-res[,10])), mean(-(true-res[,11])), mean(-(true-res[,12])))
        sb_maf_pop <- c(mean(-(pop-res[,10])), mean(-(pop-res[,11])), mean(-(pop-res[,12])))

	cat(S, D, nsites, "RMSE", "sample", "MLE_(unknown)", rmse_maf[1], "\n")
	cat(S, D, nsites, "RMSE", "sample", "MLE_(known)", rmse_maf[2], "\n")
	cat(S, D, nsites, "RMSE", "sample", "Expectation_(known)", rmse_maf[3], "\n")
	cat(S, D, nsites, "RMSE", "population", "MLE_(unknown)", rmse_maf_pop[1], "\n")
        cat(S, D, nsites, "RMSE", "population", "MLE_(known)", rmse_maf_pop[2], "\n")
        cat(S, D, nsites, "RMSE", "population", "Expectation_(known)", rmse_maf_pop[3], "\n")
 	cat(S, D, nsites, "bias", "sample", "MLE_(unknown)", sb_maf[1], "\n")
        cat(S, D, nsites, "bias", "sample", "MLE_(known)", sb_maf[2], "\n")
        cat(S, D, nsites, "bias", "sample", "Expectation_(known)", sb_maf[3], "\n")
        cat(S, D, nsites, "bias", "population", "MLE_(unknown)", sb_maf_pop[1], "\n")
        cat(S, D, nsites, "bias", "population", "MLE_(known)", sb_maf_pop[2], "\n")
        cat(S, D, nsites, "bias", "population", "Expectation_(known)", sb_maf_pop[3], "\n")

	}
}




