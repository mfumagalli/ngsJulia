
fin <- commandArgs(T)[1]

res <- readLines(fin)

mle <- as.numeric(strsplit(res[5], split=":")[[1]][2])

tmp <- strsplit(strsplit(res[3], split="\\[")[[1]][2],split="\\]")[[1]]

likes <- unlist(strsplit(tmp, split="; "))

lrts <- c()

for (i in 1:length(likes)) lrts <- c(lrts, diff(sort(as.numeric(strsplit(likes[i], split=" ")[[1]]), dec=T)[2:1])*2)

cat("LRT-ploidy:\t", min(lrts), "\n")

lrt_ane <- min(as.numeric(strsplit(unlist(strsplit(strsplit(res[6], split="\\[")[[1]][2],split="\\]")), split=", ")[[1]]))

cat("LRT-aneuploidy:\t", lrt_ane, "\n")




