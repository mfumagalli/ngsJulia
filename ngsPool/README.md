# ngsPool

Estimation of allele frequency (and more) from pooled sequencing data using `ngsJulia`.
`ngsPool` implements several estimators of allele frequencies from pooled NGS data.
It provides additional scripts for estimators of the site frequency spectrum (SFS) and for association tests.

To showcase its use, this tutorial will simulate some pooled sequencing data and demonstrate the various options and possible analyses implemented in `ngsPool`.
Before starting, please specify paths to both Julia language and ngsJulia (yours could be different):
```
JULIA=~/Software/julia-1.6.1/bin/julia
NGSJULIA=~/Software/ngsJulia
```

## Simulate pooled NGS data 

We can simulate NGS data from a pooled sequencing experiment using a script provided in `ngsJulia`.
We can explore its options:
```
Rscript $NGSJULIA/simulMpileup.R --help
```

Let's assume that we wish to simulate 10 diploid genomes, 1000 base pairs each with an average sequencing depth of 20 and base quality of 20 in Phred score. Samples come from a population of 10,000 effective size under constant-size evolution.
We can do that by running:
```
Rscript $NGSJULIA/simulMpileup.R --out test.txt --copy 2x10 --sites 1000 --depth 20 --qual 20 --ksfs 1 --ne 10000 --pool | gzip > test.mpileup.gz
```
and explore the files generated with:
```
ls test.*
```
The file `test.txt` contains the true genotypes while `test.mpileup.gz` is a gzipped [mpileup](http://www.htslib.org/doc/samtools-mpileup.html) file containing information on sequencing data 
	
## Estimate allele frequencies with SNP calling and unknown sample size

Let's explore of `ngsPool` can be used to estimate allele frequencies.
We can retrieve a list of all available options by typing:
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

