
options(scipen=999)

library("getopt")

# http://www.inside-r.org/packages/cran/getopt/docs/getopt.package
spec=matrix(c(
	'ksfs', 'k', 2, "double", "coeff. for shape of SFS default [1.0]",
	'ne', 'n', 2, "integer", "effective population size [default 10,000]",
	'snpcall', 's', 0, "logical", "flag id snps were called",
	'panc', 'p', 2, "double", "prob of ancestral being correct, if 0.5 is folded [default], if<0 compute",
	'help', 'h', 0, "logical", "print help message"
	), byrow=TRUE, ncol=5)
opt = getopt(spec)

# help
if ( !is.null(opt$help) ) {
	write.table(spec, sep="\t", quote=F, col.names=F, row.names=F)
	q(status=1)
}

if (is.null(opt$ksfs)) opt$ksfs=1.0
if (is.null(opt$ne)) opt$ne=1e4
if (is.null(opt$panc)) opt$panc=0.5
if (is.null(opt$snpcall)) opt$snpcall=FALSE

args=commandArgs(T)
Ne=opt$ne
K=opt$ksfs
snpcall=opt$snpcall
panc=opt$panc
rm(opt)

# all polymoprhic in the population!
ee=(1/(1:(Ne-1))^(1/K)); ee=ee/sum(ee);
pder=c(0, ee, 0);

# this is the expected derived allele frequency
q=weighted.mean(seq(0,Ne,1), pder)/Ne
p=1-q

# if panc=0.5 then it is folded
if (panc<0) panc=sum(pder[1:(floor(Ne/2)+1)])

if (K>0) {

	if (snpcall==0) {

		cat(c(panc, 1-panc), sep="\t")
		cat("\n")

		cat(c(p,q), sep="\t")
		cat("\n")

		cat(c(p^2,2*p*q,q^2), sep="\t")
		cat("\n")

		cat(c(p^3,3*p^2*q,3*p*q^2,q^3), sep="\t")
		cat("\n")

		cat(c(p^4,4*p^3*q,6*p^2*q^2,4*p*q^3,q^4), sep="\t")
		cat("\n")

		cat(c(p^5,5*p^4*q,10*p^3*q^2,10*p^2*q^3,5*p*q^4,q^5), sep="\t")
		cat("\n")

		cat(c(p^6, 6*p^5*q, 15*p^4*q^2, 20*p^3*q^3, 15*p^2*q^4, 6*p*q^5, q^6 ), sep="\t")
		cat("\n")

		freq <- 1-p
		cat( c( (1-freq)^7, 7*(1-freq)^6*freq, 21*(1-freq)^5*freq^2, 35*(1-freq)^4*freq^3, 35*(1-freq)^3*freq^4, 21*(1-freq)^2*freq^5, 7*(1-freq)*freq^6, freq^7 ), sep="\t")
		cat("\n")

       		cat( c( (1-freq)^8, 8*(1-freq)^7*freq, 28*(1-freq)^6*freq^2, 56*(1-freq)^5*freq^3, 70*(1-freq)^4*freq^4, 56*(1-freq)^3*freq^5, 28*(1-freq)^2*freq^6, 8*(1-freq)*freq^7, freq^8 ), sep="\t")
                cat("\n")

	} else {

		cat(c(panc, 1-panc), sep="\t")
		cat("\n")

		cat(c(0,0), sep="\t")
		cat("\n")

		cat(c(0,2*p*q,0), sep="\t")
		cat("\n")

		cat(c(0,3*p^2*q,3*p*q^2,0), sep="\t")
		cat("\n")

		cat(c(0,4*p^3*q,6*p^2*q^2,4*p*q^3,0), sep="\t")
		cat("\n")

		cat(c(0,5*p^4*q,10*p^3*q^2,10*p^2*q^3,5*p*q^4,0), sep="\t")
		cat("\n")

		cat(c(0, 6*p^5*q, 15*p^4*q^2, 20*p^3*q^3, 15*p^2*q^4, 6*p*q^5, 0), sep="\t")
		cat("\n")

		freq <- 1-p
                cat( c( 0, 7*(1-freq)^6*freq, 21*(1-freq)^5*freq^2, 35*(1-freq)^4*freq^3, 35*(1-freq)^3*freq^4, 21*(1-freq)^2*freq^5, 7*(1-freq)*freq^6, 0 ), sep="\t")
                cat("\n")

                cat( c( 0, 8*(1-freq)^7*freq, 28*(1-freq)^6*freq^2, 56*(1-freq)^5*freq^3, 70*(1-freq)^4*freq^4, 56*(1-freq)^3*freq^5, 28*(1-freq)^2*freq^6, 8*(1-freq)*freq^7, 0 ), sep="\t")
                cat("\n")

	}

} else { # K=0

	if (snpcall==0) {

		cat(c(1,1), sep="\t")
		cat("\n")

		cat(c(1,1), sep="\t")
		cat("\n")

		cat(c(1,1,1), sep="\t")
		cat("\n")

		cat(c(1,1,1,1), sep="\t")
		cat("\n")

		cat(c(1,1,1,1,1), sep="\t")
		cat("\n")

		cat(c(1,1,1,1,1,1), sep="\t")
		cat("\n")

		cat(c(1,1,1,1,1,1,1), sep="\t")
		cat("\n")

		cat(rep(1,7), sep="\t")
		cat("\n")

		cat(rep(1,8), sep="\t")
                cat("\n")

	} else {

		cat(c(1,1), sep="\t")
		cat("\n")

		cat(c(0,0), sep="\t")
		cat("\n")

		cat(c(0,1,0), sep="\t")
		cat("\n")

		cat(c(0,1,1,0), sep="\t")
		cat("\n")

		cat(c(0,1,1,1,0), sep="\t")
		cat("\n")

		cat(c(0,1,1,1,1,0), sep="\t")
		cat("\n")

		cat(c(0,1,1,1,1,1,0), sep="\t")
		cat("\n")

		cat(c(0,1,1,1,1,1,1,0), sep="\t")
                cat("\n")

		cat(c(0,1,1,1,1,1,1,1,0), sep="\t")
                cat("\n")

	}


}
