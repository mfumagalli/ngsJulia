
args <- commandArgs(T)
fin <- args[1]
rm(args)

res <- read.table(file=fin, sep="\t", head=T, stringsAsFactors=F)

ind <- which.max(res$likeRatio)

lams <- (res$mDepth[ind]/res$ploidy[ind])*seq(1,6)

for (i in 1:nrow(res)) {

	priors <- rep(NA, 6)

	priors <- dpois(x=round(res$mDepth[i]), lambda=lams)
	priors <- priors/sum(priors)

	cat(priors, sep="\t")
	cat("\n")
}



