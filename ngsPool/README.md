# ngsPool
Estimation of allele frequency from pooled-sequencing data

## Installation

        cd ngsPool
        ln -s ../ngsJulia/simulMpileup.R simulMpileup.R
        ln -s ../ngsJulia/generics.jl generics.jl
	ln -s ../ngsJulia/templates.jl templates.jl

## Simulate pool data 

```
JULIA=~/Software/julia-1.6.1/bin/julia
NGSJULIA=~/Software/ngsJulia

Rscript $NGSJULIA/simulMpileup.R --help

# e.g. simulate 10 diploids, 1000 base pairs with an average depth of 20 and base quality of 20 in phred schore, from a population of 10,000 effective size under constant-size evolution:

Rscript $NGSJULIA/simulMpileup.R --out test.txt --copy 2x10 --sites 1000 --depth 20 --qual 20 --ksfs 1 --ne 10000 --pool | gzip > test.mpileup.gz

ls test.*
```
	
## Estimate allele frequencies with SNP calling and frequentist approach

```
$JULIA $NGSJULIA/ngsPool/ngsPool.jl --help
```

```
$JULIA $NGSJULIA/ngsPool/ngsPool.jl --fin test.mpileup.gz --fout test.out.gz --lrtSnp 7.82

less -S test.out.gz
```

## Estimate allele frequency likelihoods without SNP calling

`nChroms` is equal to (ploidy x individuals), 2x8 + 3x4 = 28. If this value is set it enables the calculation of sample allele frequency likelihoods (`saf.gz` file).

```
$JULIA $NGSJULIA/ngsPool/ngsPool.jl --fin test.mpileup.gz --fout test.out.gz --nChroms 28 --fsaf test.saf.gz

less -S test.out.gz

less -S test.saf.gz
```

## Site frequency spectrum

```
Rscript $NGSJULIA/ngsPool/poolSFS.R test.saf.gz > sfs.txt

less -S sfs.txt
```

# Association test





## Assessment of its performance
Compared with Popoolation2, Snape and VarScan (https://github.com/Amend-1634/ngsPool_assessment)

