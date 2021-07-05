# ngsPloidy
Inference of ploidy from short-read sequencing data.

In this application, we show how `ngsJulia` can be used to implement a statistical test of ploidy level assignment from sequencing data.
`ngsPloidy` takes in input gzipped mpileup files and returns several metrics of per-sample ploidy assignment. It also tests for aneuploidy within the sample. 
'ngsPloidy` supports data filtering based on quality and depth values.


## Model

With the following notation:
* P: probability
* D: observed sequencing data
* O: ploidy
* A: ancestral state
`ngsPloidy` calculates the following log-likelihood for the whole sample:
```
log(P(D, G, A|O_1=y_1, O_2=y_2, ..., O_n=y_n)) = sum_{samples} sum_{sites} P(D|G,A) P(G|O,A) P(A)
```
where:
* P(D|G,A): genotype likelihood
* P(G|O,A): genotype probability (derived from the expected population allele frequency from the site frequency spectrum)
* P(A): probability of correct assignment of ancestral state
The model assumes that we observed at most two alleles. It can also calculate the per-sample log-likelihood by omitting ```sum_{samples}```.

`ngsPloidy` will output the most likely array of marginal ploidies, as well as the likelihood of all samples having the same ploidy.

In the latter assumption, we can propose a Bayesian formulation:
```
P(O_1=y, O_2=y, ..., O_n=y|D) = P(D|O) P(O) / P(D)
```
where P(D|O) is calculated as aforementioned.

## Tutorial

### Initialisation

Let's initialise paths to Julia language and `ngsJulia`.
```
JULIA=~/Software/julia-1.6.1/bin/julia
NGSJULIA=~/Software/ngsJulia
```

First, we need a create a file containing probabilities of genotypes and of the major allele being ancestral. These probability files can be generated using the following R script:
```
Rscript $NGSJULIA/ngsPloidy/writePars.R --help
```
where:
* '-k' denotes the shape of the site frequency spectrum:
        - k=1 : constant population size
        - k>1 : population bottleneck
        - k<1 : population growth
* '-n' is the effective population size
* '-s' if a flag and specifies that SNPs have beeen called; this makes sense only if only one samples is analysed with called SNPs
* '-p' specifies how to define the probability that the ancestral state is the major allele, if 0.5 this means you assume folded data, if -1 it will compute it using the site frequency spectrum
* '-h' prints a help message.

Therefore, we can create our file as:
```
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 1 > test.pars
# without ancestral allele (probability of being correct of 0.90)
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 0.90 > test.unk.pars
# with automatic calculation of probability of misassignment of ancestral allele
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p -1 > test.auto.pars
# with folded spectrum
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 0.5 > test.fold.pars
# with SNPs only
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 0.9 -n 100000 -p 1 -s > test.snp.pars
```

As an example,
```
cat test.pars
```
outputs in the first line the probabiolity of major allele being ancestral and derived, while in the remaining lines it gives the prior genotype probabilities for diallelic sites up to a maximum ploidy of 8.


### Simulations

We can simulate sequencing data in mpileup format by specifying the ploidy of each individual and other parameters of the sequencing experiment and species.
We can use the script provided:
```
Rscript $NGSJULIA/simulMpileup.R --help
```

### Case studies

Let's investigate several case studies to appreciate how `ngsPloidy` works.
We can run
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --help
```
to see all options available.

#### Case A: 2 haploids, 2 diploids, 2 triploids, 2 tetraploids, 2 pentaploids, sequenced a 20X at the haploid level

We simulate 100 sites all polymorphic in the population.
```
Rscript $NGSJULIA/simulMpileup.R --out test.A.txt --copy 1x2,2x2,3x2,4x2,5x2 --sites 100 --depth 20 | gzip > test.A.mpileup.gz

# true data
less -S test.A.txt

# observed sequencing data
less -S test.A.mpileup.gz
```

If we assume that we know the ancestral state and it's equivalent to the reference state, we can run:
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --fpars test.pars --nSamples 10 --keepRef 1
```
The option `--keepRef` forces the reference allele to be one of the two considered alleles and it is mandatory with `--fpars``.

If `--fout` is given, the program will print some statistics for each site, including the estimated allele frequency.
```
# automatic set of probability of major allele being the ancestral state
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --fpars test.auto.pars --fout test.A.out.gz --nSamples 10 --keepRef 1

less -S test.A.out.gz
```

### Case B: 10 triploids and no output file, simulating an error rate of assigning the ancestral state of 0.10

	Rscript simulMpileup.R --out test.B.txt --copy 3x10 --sites 5000 --depth 100 --qual 20 --ksfs 1 --ne 10000 --panc 0.10 | gzip > test.B.mpileup.gz

	julia ngsPoly.jl --fin test.B.mpileup.gz --fpars test.unk.pars --nSamples 10 --thSnp -1 --ploidy 1-5

### Case C: 1 diploid, 8 triploids, 1 tetraploid with ploidy prior and folded data

	Rscript simulMpileup.R --out test.C.txt --copy 2x1,3x8,4x1 --sites 5000 --depth 100 --qual 20 --ksfs 1 --ne 10000 | gzip > test.C.mpileup.gz

	julia ngsPoly.jl --fin test.C.mpileup.gz --fpars test.fold.pars --fout testC.out.gz --nSamples 10 --thSnp -1 --ploidy 1-5 --prior 0,0.1,0.8,0.1,0

### Case D: 1 tetraploid with Ne=1e6 and experiencing population growth with SNP calling

	Rscript simulMpileup.R --out test.D.txt --copy 4x1 --sites 5000 --depth 100 --qual 20 --ksfs 0.90 --ne 100000 | gzip > test.D.mpileup.gz

	julia ngsPoly.jl --fin test.D.mpileup.gz --fpars test.snp.pars --fout test.D.out.gz --nSamples 1 --thSnp 7.82 --ploidy 1-5

Please note that when `test.snp.pars` is used haploid likelihoods cannot be calculated.
Also note that by default no SNP calling is performed.

## Real data

In practise, one should also exclude trialleic sites (or in general with more than two alleles).
This can be achieved by setting a threshold on `--thTria` (default is Inf).




