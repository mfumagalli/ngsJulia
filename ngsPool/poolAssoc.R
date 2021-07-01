
args <- commandArgs(T)

saf_cases <- read.table(args[1])
saf_controls <- read.table(args[2])

p0 <- apply(FUN=which.max, X=safs[[2]], MAR=1) # in practice should be -1, these are indexs

                for (j in 1:length(p0)) {

                        H0 <- safs[[1]][j,p0[j]]
                        HA <- max(safs[[1]][j])
                        LRT <- -2*(H0-HA)

                        cat(D,F,LRT,"\n")





