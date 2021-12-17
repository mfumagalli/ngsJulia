
# calculate confusion matrix of SNP calling

cat("S","D","F","nsites","True","Predicted","Value","F1","\n")

res <- read.table("results.txt", head=T)

for (S in sort(unique(res$S))) {

	for (D in sort(unique(res$D))) {

		sub <- res[which(res$S==S & res$D==D),]

		for (F in sort(unique(sub$F))) {		

			# SNP calling (positive is SNP)

			tn <- which(sub$F==0 & sub$est==0)
			tp <- which(sub$F==F & sub$est>0)
			fn <- which(sub$F==F & sub$est==0)
			fp <- which(sub$F==0 & sub$est>0)

			f1 <- length(tp) / (length(tp) + 0.5*(length(fp)+length(fn)))

			ns1 <- length(tn)+length(fp)
			ns2 <- length(tp)+length(fn)

			if((length(tn)+length(fp))!=(length(tp)+length(fn))) break("wrong sum check!")

			cat(S,D,F,ns1,0,0,length(tn),f1,"\n")
			cat(S,D,F,ns1,1,1,length(tp),f1,"\n")
			cat(S,D,F,ns1,1,0,length(fn),f1,"\n")
			cat(S,D,F,ns1,0,1,length(fp),f1,"\n")
		}

	}

}










