
# estimate SFS from true sample frequencies and .saf files

# fit alfa and beta of a Beta distribution

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


getExpSFS <- function(S,K) {
	ee <- (1/(1:((S*2)-1))^(1/K))
	ee/sum(ee)
}


getProbBeta2 <- function(n, alfa, beta, minfreq=0.005) {
        x <- seq(0, (n*2), 0.5 ) / (n*2)
        cum <- pbeta(x, alfa, beta)
        probs <- diff(cum[seq(2,length(x),2)])
	first <- cum[2] #- pbeta(minfreq, alfa, beta)
	last <- 1-cum[length(cum)-1] #pbeta(1-minfreq, alfa, beta)-cum[length(cum)-1]
	probs <- c(first, probs, last)
        probs / sum(probs)
}


getProbBeta <- function(n, alfa, beta) {
	x <- seq(0, (n*2), 0.5 ) / (n*2)
	cum <- pbeta(x, alfa, beta)
	probs <- diff(cum[seq(2,length(x),2)])
	probs / sum(probs)
}

getDensBeta <- function(n, alfa, beta) {
        x <- seq(1, (n*2)-1) / (n*2)
        dens <- dbeta(x, alfa, beta)
        dens / sum(dens)
}

getDensGamma <- function(n, shape, scale) {
        x <- seq(1, ((n*2)-1), 1) / (n*2)
        dens <- dgamma(x, shape=shape, scale=scale)
        dens / sum(dens)
}

getCumBeta <- function(n, alfa, beta) {

	x <- seq(1, (n*2)-1) / (n*2)
	probs <- pbeta(x, alfa, beta)
	probs
}


getFit <- function(vec,Q=NA) {
	P <- vec; #P <- P/sum(P) + th
        KL <- sum(P*log(P/Q), na.rm=T)
        KL2 <-sum(Q*log(Q/P), na.rm=T)
        MSE <- sum((Q-P)^2,na.rm=T)/length(P)
        RMSE <- sum(sqrt((Q-P)^2),na.rm=T)/length(P)
	res <- c(KL, KL2, MSE, RMSE)
	#print(res)
	res
}

calcFit <- function(par, data, metric="KL", func=getProbBeta) {

	if(is.null(dim(data))) data <- t(data)
	n <- (ncol(data)+1)/2

	Q <- func(n, par[1], par[2]) #+ th 

	metrics <- t(apply(FUN=getFit, X=data, Q=Q, MAR=1))

	res <- apply(FUN=sum, X=metrics, MAR=2, na.rm=T)

	if (metric=="KL") return(res[1]); if (metric=="KL2") return(res[2])
	if (metric=="MSE") return(res[3]); if (metric=="RMSE") return(res[4])
}


norm <- function(a,th=0.01) {
	a[which(a<0.01)] <- NA
	a/sum(a, na.rm=T)
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


for (S in c(20,50)) {

	for (D in c(0.5, 1, 2, 5)) {

		# estimated frequencies
		est <- read.table(paste("../analysis_1/Results/",S,".",D,".out.gz",sep="",collapse=""))
		# saf values
        	saf <- read.table(paste("../analysis_1/Results/",S,".",D,".saf.gz",sep="",collapse=""))
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

		# 4) SFS fit from sfs_mle
		sfs_mle_snps <- sfs_mle[2:(S*2)]
		sfs_mle_snps <- sfs_mle_snps / sum(sfs_mle_snps)
		opt_mle <- optim(par=c(1,1), fn=calcFit, data=sfs_mle_snps, lower=c(1e-4,1e-4), upper=c(1e3,1e3), method="L-BFGS-B", metric="RMSE")
		sfs_opt_mle <- c(NA, getProbBeta(S, opt_mle$par[1], opt_mle$par[2])*nsnps, NA)

		# 5) SFS fit from saf_max
                opt_max1 <- optim(par=opt_mle$par, fn=calcFit, data=saf_max1, lower=c(1e-4,1e-4), upper=c(1e3,1e3), method="L-BFGS-B", metric="RMSE")
                sfs_opt_max1 <- c(NA, getProbBeta(S, opt_max1$par[1], opt_max1$par[2])*nrow(saf_max1), NA)

		# 6) SFS fit from saf_max3
                opt_max3 <- optim(par=opt_mle$par, fn=calcFit, data=saf_max3, lower=c(1e-4,1e-4), upper=c(1e3,1e3), method="L-BFGS-B", metric="RMSE")
                sfs_opt_max3 <- c(NA, getProbBeta(S, opt_max3$par[1], opt_max3$par[2])*nrow(saf_max3), NA)

		# results
		res <- rbind(sfs_pop, sfs_sample, sfs_mle, sfs_opt_mle, sfs_opt_max1, sfs_opt_max3)

		# plot
		#barplot(res, beside=T)
		for (i in 1:nrow(res)) {
			cat(S, D, rownames(res)[i], sum(abs(res[i,]-res[1,]),na.rm=T), sum(abs(res[i,]-res[2,]), na.rm=T), res[i,], "\n")
		}

	}
}


# STILL ISSSUES WITH OPTIMIZATION, CHECK AGAIN, EVEN with rmse, do some tests what happens with 0s and NAs, max1 should work!






