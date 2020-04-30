# ngsPool
Population genetics from pool sequencing data

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

## Estimate allele frequencies (MLE only with SNP calling)

	julia ngsPool.jl --help

	julia ngsPool.jl --fin test.mpileup.gz --fout test.out.gz --thSnp 7.82

	less -S test.out.gz

## Estimate allele frequency likelihoods (without SNP calling)

	julia ngsPool.jl --fin test.mpileup.gz --fout test.out.gz --thSnp -Inf --nChroms 220 --fsaf test.saf.gz

	less -S test.out.gz
	less -S test.saf.gz


