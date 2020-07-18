# ngsPool
Estimation of minor allele frequency from pooled-sequenced Next-generation-sequencing data, which is updated and commented version.

## Installation

	git clone https://github.com/mfumagalli/ngsPool.git
        git clone https://github.com/mfumagalli/ngsJulia.git

        cd ngsPool
        ln -s ../ngsJulia/simulMpileup.R simulMpileup.R
        ln -s ../ngsJulia/generics.jl generics.jl
	ln -s ../ngsJulia/templates.jl templates.jl

## Simulate pool data

	Rscript simulMpileup.R --help

	Rscript simulMpileup.R --out test.txt --copy 2x80,3x20 --sites 1000 --depth 100 --qual 20 --ksfs 1 --ne 10000 --pool | gzip > test.mpileup.gz

	ls test.*
	
### Simulating pool data for case
	time Rscript simulMpileup_case.R --rr 2 --prevalence 0.1 --out_qqVector data2/case_control/qqVector-1-1.txt --out data2/case_control/test-case-1-1.txt --copy 2x80 --sites 1000 --depth 5 --pool | gzip > data2/case_control/mpileup-case-1-1.gz

### Simulating pool data for control
	time Rscript simulMpileup_control.R --rr 2 --prevalence 0.1 --qqVector qqVector-1-1.txt --out test-control-1-1.txt --copy 2x80 --sites 1000 --depth 5 --pool | gzip > mpileup-control-1-1.gz

     
    
##    
*input -nSamp number can generate meaningful --saf(Site frequency spectrum) output and minor allele frequency estimate based on SFS (freqMax (maximum likelihood) and freqE(expected))
## Estimate allele frequencies (MLE only with SNP calling)
Default likelihood ratio test statistic cutoff is 7.82, filtering out polymorphic sites.
	julia ngsPool.jl --help

	julia ngsPool.jl --fin test.mpileup.gz --fout test.out.gz --nSamp 220 --fsaf test.saf.gz --lrtSnp 7.82

	less -S test.out.gz

## Estimate allele frequency likelihoods (without SNP calling)
all sites.
	
	julia ngsPool.jl --fin test.mpileup.gz --fout test.out.gz --nSamp 220 --fsaf test.saf.gz --lrtSnp -Inf

	less -S test.out.gz
	less -S test.saf.gz


