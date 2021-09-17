
fin <- commandArgs(T)[1]

res <- readLines(fin)

mle <- as.numeric(strsplit(res[5], split=":")[[1]][2])

tmp <- strsplit(strsplit(res[3], split="\\[")[[1]][2],split="\\]")[[1]]

likes <- unlist(strsplit(tmp, split="; "))

nr_samples <- length(likes)

lrts <- c()

for (i in 1:length(likes)) lrts <- c(lrts, diff(sort(as.numeric(strsplit(likes[i], split=" ")[[1]]), dec=T)[2:1])*2)


inferred <- strsplit(res[4], split=":")[[1]][2]


cat("The inferred vector of ploidy levels is", inferred, "with a LR value of", min(lrts), "between the most and second most likely vectors.", "\n")
#cat("LRT-ploidy:\t", min(lrts), "\n")

lrt_ane <- min(as.numeric(strsplit(unlist(strsplit(strsplit(res[6], split="\\[")[[1]][2],split="\\]")), split=", ")[[1]]))

pv <- pchisq(lrt_ane, df=(nr_samples-1), lower.tail=F)

cat("The statistical support for rejecting the hypothesis of equal ploidy across all samples is", lrt_ane, "in LRT value with a p-value of", pv,".\n")




