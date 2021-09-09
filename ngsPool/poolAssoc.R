
args <- commandArgs(T)

saf_cases <- read.table(args[1])
saf_controls <- read.table(args[2])

p0 <- apply(FUN=which.max, X=saf_controls, MAR=1) # the MLE of frequency is p0-1, p0 are indexes

cat("##LRT","\t","log(p-value)","\n")

for (j in 1:length(p0)) {

	H0 <- saf_cases[j,p0[j]]
        HA <- max(saf_cases[j,])
        LRT <- -2*(H0-HA)
	pv <- pchisq(LRT, df=1,lower.tail=T, log.p=T)

        cat(LRT, "\t", pv, "\n")

}



