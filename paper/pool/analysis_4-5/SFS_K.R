
# estimate SFS from true sample frequencies and .saf files

calcSFS <- function(n, daf, norm=TRUE) {
	# from relative frequencies [0,1], both population and sample
	# including non-snps
	x <- seq(0,1,1/(2*n))
	sfs <- rep(0, length(x))
	# daf is [0,1]
	for (i in daf) {
		ind <- which.min(abs(i-x))
		sfs[ind] <- sfs[ind] + 1
	}
	if (norm) sfs <- sfs/sum(sfs)
	sfs
}

getExpSFS <- function(K,S=0) {
	ee <- (1/(1:((S*2)-1))^(1/K))
	ee/sum(ee)
}

getFit <- function(vec,Q=NA,th=1e-8) {
	P <- vec + th; #P <- P/sum(P) + th
	Q <- Q + th
        KL <- sum(P*log(P/Q), na.rm=T)
        KL2 <-sum(Q*log(Q/P), na.rm=T)
        MSE <- sum((Q-P)^2,na.rm=T)/length(P)
        RMSE <- sum(sqrt((Q-P)^2),na.rm=T)/length(P)
	res <- c(KL, KL2, MSE, RMSE)
	res
}

calcFit <- function(par, data, metric="KL", func=getExpSFS, th=1e-8) {

	if(is.null(dim(data))) data <- t(data)
	n <- (ncol(data)+1)/2

	Q <- func(S=n, K=par[1]) 

	metrics <- t(apply(FUN=getFit, X=data, Q=Q, MAR=1, th=th))

	res <- apply(FUN=sum, X=metrics, MAR=2, na.rm=T)

	if (metric=="KL") return(res[1]); if (metric=="KL2") return(res[2])
	if (metric=="MSE") return(res[3]); if (metric=="RMSE") return(res[4])
}


norm <- function(a,th=0.05) {
	a[which(a<th)] <- 0
	a/sum(a)
}

maxim <- function(a) {
	b <- rep(0, length(a))
	b[which.max(a)] <- 1
	b
}

maxim3 <- function(a) {
        b <- rep(0, length(a))
	ind <- which.max(a)
	inds <- c(ind-1, ind, ind+1)
	inds <- inds[which(inds>0 & inds<=length(a))]
        b[inds] <- a[inds]
        b/sum(b)
}

fout <- "results.4.txt"
cat("S", "D", "estimation", "reference", "rmse", "\n", file=fout)

foutK <- "results.4.K.txt"
cat("S","D","estimation","value", "\n", file=foutK)


for (S in c(20,50,100)) {

	fsfs <- paste("results.4.",S,".sfs.txt",sep="",collapse="")
	cat("", file=fsfs)

	for (D in c(1, 2, 5)) {

		cat(S,D,"\n",file=stderr())

		# estimated frequencies
		est <- read.table(paste("../analysis_1/Results/",S,".",D,".out.gz",sep="",collapse=""))
		# saf values
        	saf <- exp(read.table(paste("../analysis_1/Results/",S,".",D,".saf.gz",sep="",collapse="")))
		# simulated true frequencies
        	sim <- read.table(paste("../analysis_1/Data/",S,".",D,".txt",sep="",collapse=""))

		# take only values with data, so rows in estimated
                sim <- sim[est[,2],]
                nsites <- nrow(sim)

		# 0) SFS populationfrom theory
                #sfs_theory <- c(NA,(1/(1:((S*2)-1))),NA)
                #sfs_theory <- sfs_theory / sum(sfs_theory)

		# 1) SFS population
		daf_pop <- sim[,5]
		sfs_pop <- calcSFS(S, daf_pop, F)
		nsnps <- sum(sfs_pop[2:(S*2)])

		# 2) SFS sample
        	daf_sample <- sim[,(S+6)]
		sfs_sample <- calcSFS(S, daf_sample, F)

		# 3) SFS mle (max saf) 
		daf_mle <- (apply(X=saf, MAR=1, FUN=which.max)-1) / (S*2)
		sfs_mle <- calcSFS(S, daf_mle, F)

		## only SNPs
		ind <- which(daf_mle>0 & daf_mle<1)
		saf_max1 <- t(apply(X=saf[ind,2:(S*2)], FUN=maxim, MAR=1))
		saf_max3 <- t(apply(X=saf[ind,2:(S*2)], FUN=maxim3, MAR=1))
		saf_th <- t(apply(X=saf[ind,2:(S*2)], FUN=norm, MAR=1, th=0.05))
		nsnps_sample <- length(ind)

		# record Ks
		Ks <- c()

		# 4) SFS fit from sfs_mle
		sfs_mle_snps <- sfs_mle[2:(S*2)]
		sfs_mle_snps <- sfs_mle_snps / sum(sfs_mle_snps)
		opt_mle <- optim(par=c(1), fn=calcFit, data=sfs_mle_snps, lower=c(0), upper=c(1e2), method="L-BFGS-B", metric="KL")
		sfs_opt_mle <- c(NA, getExpSFS(S=S, K=opt_mle$par[1])*nsnps_sample, NA)
		Ks <- c(Ks, opt_mle$par[1])

		# 5) SFS fit from saf_max
                opt_max1 <- optim(par=opt_mle$par, fn=calcFit, data=saf_max1, lower=c(0), upper=c(1e2), method="L-BFGS-B", metric="KL")
                sfs_opt_max1 <- c(NA, getExpSFS(S=S, K=opt_max1$par[1])*nsnps_sample, NA)
		Ks <- c(Ks, opt_max1$par[1])

		# 6) SFS fit from saf_max3
                opt_max3 <- optim(par=opt_mle$par, fn=calcFit, data=saf_max3, lower=c(0), upper=c(1e2), method="L-BFGS-B", metric="KL")
                sfs_opt_max3 <- c(NA, getExpSFS(S=S, K=opt_max3$par[1])*nsnps_sample, NA)
		Ks <- c(Ks, opt_max3$par[1])

		# 7) SFS fit from saf_th
                opt_th <- optim(par=opt_mle$par, fn=calcFit, data=saf_th, lower=c(0), upper=c(1e2), method="L-BFGS-B", metric="KL", th=1e-6)
                sfs_opt_th <- c(NA, getExpSFS(S=S, K=opt_th$par[1])*nsnps_sample, NA)
		Ks <- c(Ks, opt_th$par[1])

		# results
		res <- rbind(sfs_pop, sfs_sample, sfs_mle, sfs_opt_mle, sfs_opt_max1, sfs_opt_max3, sfs_opt_th)

		resnames <- c("pop", "sample", "count", "fit_count", "fit_afl_max1", "fit_afl_max3", "fit_afl")

		Knames <- c("fit_count", "fit_afl_max1", "fit_afl_max3", "fit_afl")

		# plot
		#barplot(res, beside=T)
		for (i in 3:nrow(res)) {
			cat(S, D, resnames[i], "pop", sum(abs(res[i,2:(S*2)]-res[1,2:(S*2)]),na.rm=T), "\n", file=fout, append=T)
			cat(S, D, resnames[i], "sample", sum(abs(res[i,2:(S*2)]-res[2,2:(S*2)]), na.rm=T), "\n", file=fout, append=T)
		}
		for (i in 1:nrow(res)) {
			cat(res[i,], "\n", file=fsfs, append=T)
		}
		# Ks
		for (i in 1:length(Ks)) {
                        cat(S, D, Knames[i], Ks[i], "\n", file=foutK, append=T)
                }


	}
}







