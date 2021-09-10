# ngsPloidy

Inference of ploidy from short-read NGS data using `ngsJulia.`
`ngsPloidy` takes gzipped mpileup files as input and returns several metrics of per-sample ploidy assignment. 
It calculates the likelihood of NGS given each tested ploidy, with genotype probabilities calculated from their likelihoods and expected (or estimated) allele frequencies.
`ngsPloidy` outputs the most likely vector of marginal ploidies, as well as the likelihood of all samples having the same ploidy.
The package also tests for aneuploidy within the sample, i.e. it provides statistical support for the ploidy being different among all tested samples.
`ngsPloidy` supports data filtering based on quality and depth values.

## Initialisation

Let's initialise paths to Julia language and `ngsJulia` (yours could be different):
```
JULIA=~/Software/julia-1.6.1/bin/julia
NGSJULIA=~/Software/ngsJulia
```

First, we need a create a file containing genotype probabilities. 
We also need to provide the probability of the major allele being ancestral, as this information will be used in case of limited sample size.
These probability files can be generated using the following R script:
```
Rscript $NGSJULIA/ngsPloidy/writePars.R --help
```
where:
* '-k' denotes the shape of the site frequency spectrum:
        - k=1 : constant population size
        - k>1 : population bottleneck
        - k<1 : population growth
* '-n' is the effective population size
* '-s' if a flag and specifies that SNPs have beeen called; this is useful only if only one sample is analysed with called SNPs
* '-p' specifies how to define the probability that the ancestral state is the major allele, e.g., if 0.5 this means you assume folded data, if -1 it will compute it using the site frequency spectrum
* '-h' prints a help message.

Therefore, depending on our settings, we can create our probability files as:
* with known allelic polarisation:
```
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 1 > test.pars
```
* with uncertain assignment of ancestral alleles (probability of being correct of 0.90)
```
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 0.90 > test.unk.pars
```
* with automatic calculation of probability of misassignment of ancestral allele
```
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p -1 > test.auto.pars
```
* with folded spectrum (unknown allelic polarisation)
```
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 0.5 > test.fold.pars
```
* or with SNPs only:
```
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 0.9 -n 100000 -p 1 -s > test.snp.pars
```

As an example,
```
cat test.pars
```
outputs in the first line the probability of major allele being ancestral and derived, while in the remaining lines it shows the probabilities for diallelic genotypes for each tested ploidy (up to ploidy 8).

## Simulations

We can simulate sequencing data in mpileup format by specifying the ploidy of each individual and other parameters of the sequencing experiment and species.
We can use the script provided:
```
Rscript $NGSJULIA/simulMpileup.R --help
```

## Case studies

Let's investigate several case studies to understand how `ngsPloidy` works.
We can run
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --help
```
to see all options available.

#### Case A: 2 haploids, 2 diploids, 2 triploids, 2 tetraploids, 2 pentaploids, sequenced at 20X at haploid level

We simulate 100 sites all polymorphic in the population with:
```
Rscript $NGSJULIA/simulMpileup.R --out test.A.txt --copy 1x2,2x2,3x2,4x2,5x2 --sites 100 --depth 20 | gzip > test.A.mpileup.gz
```
The true data is contained in this file:
```
less -S test.A.txt
```
while the simulated sequencing data is accessible with:
```
# observed sequencing data
less -S test.A.mpileup.gz
```

If we assume that we know the ancestral state and it is equivalent to the reference allele, we can run:
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --fpars test.pars --nSamples 10 --keepRef 1
```
The option `--keepRef` forces the reference allele to be one of the two considered alleles and it is mandatory with `--fpars``.

Results are printed on the screen and show:
* nr of analysed sites: vector of sites that passed filtering for each sample
* log-likelihoods of per-sample ploidies: a matrix of nr_sample X nr_ploidies with the log-likelihood of each sample having a certain ploidy (rows are separated by ;)
* #MLE vector of ploidies: the vector if individual maximum likelihood estimates of ploidy for each sample
* #log-likelikehood of MLE vector of ploidies: the log-likelihood of the above vector of estimated ploidies
* #LRT of aneuploidy: difference between MLE vector of ploidies and the log-likliehood of all samples having the same ploidy, calculated for all tested ploidies

In this example, we can see that the vector of estimated ploidy levels is equal to the simulated one.
We further observe that the LRT values for aneuploidy are high, further suggesting variation in ploidy among samples.

----------------------------------------

If `--fout` is given (optional), `ngsPloidy` prints some statistics for each site, including the estimated allele frequency:
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --fpars test.pars --fout test.A.out.gz --nSamples 10 --keepRef 1
```
with the output file accessible with:
```
less -S test.A.out.gz
```
Specifically, this file reports:
* chrom: chromosome
* pos: position
* ref: reference allele
* depth: sequencing depth
* ref/anc: reference or ancestral allele (inferred)
* alt/der: alternate or derived allele (inferred)
* lrtSnp: LRT for the site being a SNP
* lrtBia: LRT for the site being biallelic
* lrtTria: LRT for the site being triallelic
* aaf: estimated alternate or ancestral allele frequency

---------------------------------------------

If `-fglikes` is given (optional), the program ouputs the per-site genotype likelihoods for each tested ploidy.
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --fpars test.pars --fglikes test.A.glikes.gz --nSamples 10 --keepRef 1
```
with the output file accessible with:
```
less -S test.A.glikes.gz
```

--------------------------

We can also change the genotype probabilities in input. For instance, we can impose an automatic setting of the probability of major allele being the ancestral state with:
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --fpars test.auto.pars --fout test.A.out.gz --nSamples 10 --keepRef 1
```

---------------------------

One additional possibility is to not impose any genotype probability based on the site frequency spectrum (with `--fpars`) and use a uniform probability distribution instead (`--unif 1`).
With this option you are required that there are is no missing data at each sample (i.e. `--minSamples` should be equal to `--nSamples`), otherwise the program will throw an error.

```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --unif 1 --nSamples 10 --minSamples 10
```

---------------------------

Finally, we can even call genotypes (`--callGeno`) and consider their likelihoods only in the ploidy estimation with:
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --unif 1 --callGeno 1 --nSamples 10 --minSamples 10
```

Please note that we obtained the same inferred vector of ploidy levels in all these different options.
However, results may vary significantly between different settings with lower sequencing depth or lower quality.
For instance, calling genotypes is not recommended for low-depth data.


### Case B: 10 triploids and no output file, simulating an error rate of assigning the ancestral state of 0.10

We can simulate this scenario with:
```
Rscript $NGSJULIA/simulMpileup.R --out test.B.txt --copy 3x10 --sites 100 --depth 20 --panc 0.1 | gzip > test.B.mpileup.gz
```

As there is uncertainty in the allelic polarisation, we can incorporate such error in the genotype probability file previously calcolated.
As a further illustration, we will filter out bases with a quality lower than 15 in Phred score.
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl	--fin test.B.mpileup.gz --fpars test.unk.pars --nSamples 10 --keepRef 1 --minQ 15
```

As you can evinced from the resulted, the test for aneuploidy is rejected as the most likely vector of ploidies supports all samples being triploid. 

## Case C: 1 diploid, 8 triploids, 1 tetraploid with and folded allele frequencies

This scenario can be simulated with:
```
Rscript $NGSJULIA/simulMpileup.R --out test.C.txt --copy 2x1,3x8,4x1 --sites 100 --depth 20 | gzip > test.C.mpileup.gz
```

As the data is folded (we have no information on which allele is ancestral or derived), we can use the appropriate genotype probability previously calculated:
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.C.mpileup.gz --fpars test.fold.pars --nSamples 10 --keepRef 1
```

From the results, aneuploidy is still statistically supported.

### Case D: 1 tetraploid with Ne=1e6, experiencing population growth with SNP calling

In this case we have only one sample and we will attempt to infer its ploidy for called SNPs.
Let's assume we have 1000 sites polymorphic in the population. As we have one sample, we assume it is sequenced at higher depth.
This scenario can be simulated with (note that, as an illustration, we are not generating an output file for the real data):
```
Rscript $NGSJULIA/simulMpileup.R --copy 4x1 --sites 1000 --depth 30 --ksfs 0.90 --ne 100000 | gzip > test.D.mpileup.gz
```

As we wish to call SNPs, we need to indicate the appropriate file for genotype probabilities and a threshold for SNP calling (`--thSnp`).
The latter is in X^2 score value where, for instance, 6.64 is equivalent to a p-value of 0.01.

```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.D.mpileup.gz --fpars test.snp.pars --keepRef 1 --nSamples 1 --thSnp 6.64
```
Please note that when `test.snp.pars` is used haploid likelihoods cannot be calculated.
Also note that in this case the number of analysed sites is less that the number of simulated sites because only SNPs are considered.

## Further options

Other options in `ngsPloidy` can be retrieved by:
```
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --help
```
Several filtering options on base quality, depth and proportion of minor reads are available. Likewise, options to include only bialleic sites (i.e. exclude triallelic and multiallelic sites) can be used.
Results may vary depending on the filtering options and users are encourage to consider how their inferences are robust to the data processing pipeline.




