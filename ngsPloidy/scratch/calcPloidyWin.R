
# sliding windows scan ploidy

args <- commandArgs(T)
fname <- args[1] # .out.gz
winSize <- as.numeric(args[2])
stepSize <- as.numeric(args[3])
thCall <- args[4] # pv
minMaf <- as.numeric(args[5])
fpriors <- args[6]

res <- read.table(fname)

starts <- seq(min(res$V2), max(res$V2), stepSize)
ends <- starts+winSize-1
lastWin <- which(ends>=max(res$V2))[1]
starts <- starts[1:lastWin]
ends <- ends[1:lastWin] 

nsites <- nsnps <- ploidy <- likeRatio <- mdepth <- sddepth <- pprobs <- c()

if (fpriors=="NULL") {
	thCall <- dchisq((1-as.numeric(thCall)),1) # pv LRT
	cat("start", "ends", "nsites", "nsnps", "mDepth", "sdDepth", "likeRatio", "ploidy", sep="\t")
} else {
	thCall <- 1-as.numeric(thCall)
	cat("start", "ends", "nsites", "nsnps", "mDepth", "sdDepth", "postProb*", "ploidy", sep="\t")
	priors <- read.table(fpriors, head=F, stringsAsFact=F, sep="\t")
}
cat("\n")

rm(args)

for (w in 1:length(starts)) {

	ind <- which(res$V2>=starts[w] & res$V2<=ends[w])
	nsites <- length(ind)
	ind <- ind[which(res$V10[ind]>=minMaf & res$V10[ind]<=(1-minMaf))]
	nsnps <- length(ind)
	mdepth <- mean(res[ind,4], na.rm=T)
	sddepth <- sd(res[ind,4], na.rm=T)

	ploidy <- NA

	if (length(ind)>0 & nsnps>0) {
	
		# test for evidence of ploidy
		likes <- as.numeric(apply(X=res[ind,11:16], FUN=sum, MAR=2))

		if (fpriors=="NULL") {
			likeRatio <- abs(diff(sort(likes, dec=T)))[1]
			if (likeRatio >= thCall) ploidy <- which.max(likes)
		} else {
			pprobs <- as.numeric(likes + log(priors[w,]))
			likeRatio <- abs(diff(sort(pprobs, dec=T)))[1]
			#pp <- pp/sum(pp)
			#if (max(pp) >= thCall) ploidy <- which.max(pp)
			ploidy <- which.max(pprobs)
		}
	}
	
	if (fpriors=="NULL") {
		cat(starts[w], ends[w], nsites, nsnps, mdepth, sddepth, likeRatio, ploidy, sep="\t")
	} else {
		cat(starts[w], ends[w], nsites, nsnps, mdepth, sddepth, likeRatio, ploidy, sep="\t")
	}
	cat("\n")
			
}


