# ngsPloidy

`nsgPloidy` implements a method to infer ploidy level from short-read NGS data using `ngsJulia.`
The model is described in the related [paper](https://f1000research.com/articles/11-126).
In a nutshell, the algorithm calculates the likelihood of ploidy levels from genotype likelihoods and genotype priors. The latter can be estimated either from the data or, with limited sample size, from a population genetic model.

--------------------------------------------

`ngsPloidy` takes gzipped mpileup files as input and returns several metrics of per-sample ploidy assignment. 
It calculates the likelihood of NGS data given each tested ploidy, with genotype probabilities calculated from their likelihoods and expected (or estimated) allele frequencies.
`ngsPloidy` outputs the most likely vector of marginal ploidies, as well as the likelihood of all samples having the same ploidy.
It also tests for multiploidy within the sample, i.e. it provides statistical support for the ploidy being different among all tested samples.
`ngsPloidy` supports data filtering based on quality and depth values.

Throughout these examples, we assume that we defined an environment variable `NGSJULIA` that points to the installation path.
Assuming that `ngsJulia` is installed in a folder name `~/Software`, this can be achieved temporarily with
```bash
NGSJULIA=~/Software/ngsJulia
```
or permanently with
```bash
sudo nano ~/.bashrc
export NGSJULIA=~/Software/ngsJulia
source ~/.bashrc
```

## Case studies

Let's investigate several case studies to understand how `ngsPloidy` works.
We can run
```bash
julia $NGSJULIA/ngsPloidy/ngsPloidy.jl --help
```
to see all options available.


### Case A: 4 haploids, 4 diploids, 4 triploids, 4 tetraploids, 4 pentaploids

We can simulate sequencing data in mpileup format by specifying the ploidy of each individual and other parameters of the sequencing experiment and species.
We can use the script provided (which requires `getopt` package):
```bash
Rscript $NGSJULIA/simulMpileup.R --help
```

We simulate 100 sites all polymorphic in the population with:
```bash
Rscript $NGSJULIA/simulMpileup.R --out test.A.txt --copy 1x4,2x4,3x4,4x4,5x4 --sites 100 --depth 20 | gzip > test.A.mpileup.gz
```
The true data is contained in this file:
```bash
zcat test.A.txt | less -S
```
where each column indicates:  
-identifier of contig (named after the --copy option)  
-position  
-reference allele (set to A)  
-alternate allele (set to C)  
-population allele frequency  
-genotype for each sample  
-sampled allele frequency  

The simulated observed sequencing data is accessible with:
```bash
zcat test.A.mpileup.gz | less -S
```
and it is formatted as standard [mpileup](http://www.htslib.org/doc/samtools-mpileup.html) file.

If we assume we have enough sample size, we can use the estimated allele frequency to calculate genotype probabilities at each site.
In this case, we can simply infer ploidy levels with:
```bash
julia $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.A.mpileup.gz --nSamples 20 > test.A.out
```
and access the output file as:
```bash
cat test.A.out
```

Results show:  
-nr of analysed sites: vector of sites that passed filtering for each sample  
-log-likelihoods of per-sample ploidies: a matrix of nr\_sample X nr\_ploidies with the log-likelihood of each sample having a certain ploidy (rows are separated by ;)  
-MLE vector of ploidies: the vector of individual maximum likelihood estimates of ploidy for each sample  
-log-likelikehood of MLE vector of ploidies: the log-likelihood of the above vector of estimated ploidies  
-LRT of multiploidy: difference between MLE vector of ploidies and the log-likelihood of all samples having the same ploidy, calculated for all tested ploidy levels  

Note that by default all ploidy levels from 1 to 8 are tested, and therefore for each sample 8 likelihoods are calculated and reported.
To extract the interpretation of these results, we can run the following script: 
```bash
Rscript $NGSJULIA/ngsPloidy/ploidyLRT.R --in test.A.out --out test.A.ploidy.out
cat test.A.ploidy.out
```
which outputs the most likely ploidy with its statistical support and the result for the test of multiploidy.
In this example, we can see that the vector of estimated ploidy levels is equal to the simulated one.
We further observe that the LRT value for multiploidy are high, further suggesting variation in ploidy among samples.
In fact, the most likely vector of equal ploidy levels is 4 but it has less support than the inferred vector of multiploidy.


### Case B: 2 haploids, 2 diploids, 2 triploids, 2 tetraploids, 2 pentaploids

This example is similar to Case A but with half of the samples.
We simulate this scenario with:
```bash
Rscript $NGSJULIA/simulMpileup.R --copy 1x2,2x2,3x2,4x2,5x2 --sites 100 --depth 20 | gzip > test.B.mpileup.gz
```

With limited sample size we can use a different estimation of population allele frequency which will be constant across all sites.
In case of limited sample size, firstly we need a create a file containing genotype probabilities.
We also need to provide the probability of the major allele being ancestral, as this information will be used in case of limited sample size.
These probability files can be generated using the following R script (which requires `getopt` package):
```bash
Rscript $NGSJULIA/ngsPloidy/writePars.R --help
```
Further documentation is available [here](https://ngsjulia.readthedocs.io/en/latest/aux/).

If we assume that we know the ancestral state and it is equivalent to the reference allele, we can run:
```bash
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 1 > test.pars

julia $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.B.mpileup.gz --fpars test.pars --nSamples 10 --keepRef 1 > test.B.out

Rscript $NGSJULIA/ngsPloidy/ploidyLRT.R --in test.B.out --out test.B.ploidy.out
```
The option `--keepRef` forces the reference allele to be one of the two considered alleles and it is mandatory with `--fpars`.

As an example,
```bash
cat test.pars
```
outputs in the first line the probability of major allele being ancestral and derived, while in the remaining lines it shows the probabilities for diallelic genotypes for each tested ploidy (up to ploidy 8).

----------------------------------------

If `--fout` is given (optional), `ngsPloidy` prints some statistics for each site, including the estimated allele frequency:
```bash
julia $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.B.mpileup.gz --fpars test.pars --fout test.B.out.gz --nSamples 10 --keepRef 1
```
with the output file accessible with:
```bash
zcat test.B.out.gz | less -S
```
Specifically, this file reports:  
-chrom: chromosome  
-pos: position  
-ref: reference allele  
-depth: sequencing depth  
-ref/anc: reference or ancestral allele (inferred)  
-alt/der: alternate or derived allele (inferred)  
-lrtSnp: LRT for the site being a SNP  
-lrtBia: LRT for the site being biallelic  
-lrtTria: LRT for the site being triallelic  
-aaf: estimated alternate or ancestral allele frequency  

Note that in this case results are printed on the screen.
Also note that if `--fpars` is not provided, `--fout` will show the estimate of the minor allele frequency (maf) instead.

--------------------------

We can also change the genotype probabilities in input. 
For instance, we can impose an automatic setting of the probability of major allele being the ancestral state with:
```bash
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p -1 > test.auto.pars

julia $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.B.mpileup.gz --fpars test.auto.pars --fout test.B.out.gz --nSamples 10 --keepRef 1
```

---------------------------

One additional possibility is to not impose any genotype probability based on the site frequency spectrum (with `--fpars`) and use a uniform probability distribution instead (`--unif 1`).
With this option you are required that there are is no missing data at each sample (i.e. `--minSamples` should be equal to `--nSamples`), otherwise the program will throw an error.

```bash
julia $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.B.mpileup.gz --unif 1 --nSamples 10 --minSamples 10
```

---------------------------

Finally, we can even infer ploidy by assigning genotypes (`--callGeno`) and consider their likelihoods only in the ploidy estimation with:
```bash
julia $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.B.mpileup.gz --unif 1 --callGeno 1 --nSamples 10 --minSamples 10
```

Please note that we obtained the same inferred vector of ploidy levels in all these different options.
However, results may vary significantly between different settings depending on sequencing depth, sample size, and number of sites.
For instance, calling genotypes is not recommended for low-depth data.


### Case C: 10 triploids and no output file, simulating an error rate of assigning the ancestral state of 0.10

We can simulate this scenario with:
```bash
Rscript $NGSJULIA/simulMpileup.R --out test.C.txt --copy 3x10 --sites 100 --depth 20 --panc 0.1 | gzip > test.C.mpileup.gz
```

As there is uncertainty in the allelic polarisation (and limited sample size), we can incorporate such error in the genotype probability file previously calcolated.
As a further illustration, we will filter out bases with a quality lower than 15 in Phred score.
```bash
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 0.90 > test.unk.pars

julia $NGSJULIA/ngsPloidy/ngsPloidy.jl	--fin test.C.mpileup.gz --fpars test.unk.pars --nSamples 10 --keepRef 1 --minQ 15 > test.C.out

Rscript $NGSJULIA/ngsPloidy/ploidyLRT.R test.C.out
```
As we can evince from the results, we fail to reject the null hypothesis of equal ploidy across all samples, as desired given the simulated scenario.

---------------------------------------------

If `-fglikes` is given (optional), the program ouputs the per-site genotype likelihoods for each tested ploidy,
```bash
julia $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.C.mpileup.gz --fpars test.unk.pars --fglikes test.C.glikes.gz --nSamples 10 --keepRef 1 --minQ 15
```
with the output file accessible with:
```bash
zcat test.C.glikes.gz | less -S
```
where genotype likelihoods (assuming diallelic variation) for all tested ploidy are provided on each line.


### Case D: 1 diploid, 8 triploids, 1 tetraploid with and folded allele frequencies

This scenario can be simulated with:
```bash
Rscript $NGSJULIA/simulMpileup.R --out test.D.txt --copy 2x1,3x8,4x1 --sites 100 --depth 20 | gzip > test.D.mpileup.gz
```

As the data is folded (we have no information on which allele is ancestral or derived), we can use the appropriate genotype probability previously calculated:
```bash
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 0.5 > test.fold.pars

julia $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.D.mpileup.gz --fpars test.fold.pars --nSamples 10 --keepRef 1 > test.D.out

Rscript $NGSJULIA/ngsPloidy/ploidyLRT.R test.D.out
```
From the results, multiploidy is statistically supported.


### Case E: 1 tetraploid with Ne=1e6, experiencing population growth with SNP calling

In this case we have only one sample and we will attempt to infer its ploidy for called SNPs.
Let's assume we have 1000 sites polymorphic in the population.
This scenario can be simulated with:
```bash
Rscript $NGSJULIA/simulMpileup.R --copy 4x1 --sites 1000 --depth 30 --ksfs 0.90 --ne 100000 | gzip > test.E.mpileup.gz
```
Note that, as an illustration, we are not generating an output file for the real data.

As we wish to call SNPs, we need to indicate the appropriate file for genotype probabilities and a threshold for SNP calling (`--lrtSnp`).
The former can be obtained with:
```bash
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 0.9 -n 100000 -p 1 -s > test.snp.pars
```
The latter is in chi-square score value where, for instance, 6.64 corresponds to a p-value of 0.01.

```bash
julia $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.E.mpileup.gz --fpars test.snp.pars --keepRef 1 --nSamples 1 --lrtSnp 6.64 > test.E.out

Rscript $NGSJULIA/ngsPloidy/ploidyLRT.R test.E.out
```
Please note that, when SNP calling in `--fpars` is used, haploid likelihoods cannot be calculated.
Also note that in this case the number of analysed sites is less that the number of simulated sites because only SNPs are considered.

## Further options

Several filtering options on base quality, depth and proportion of minor reads are available. Likewise, options to include only bialleic sites (i.e. exclude triallelic and multiallelic sites) can be used.
Results may vary depending on the filtering options and users are encourage to consider how their inferences are robust to the data processing pipeline.

All options in `ngsPloidy` can be retrieved by:
```bash
julia $NGSJULIA/ngsPloidy/ngsPloidy.jl --help

usage: ngsPloidy.jl --fin FIN [--fpars FPARS] [--fout FOUT]
                    [--fglikes FGLIKES] --nSamples NSAMPLES
                    [--ploidy PLOIDY] [--keepRef KEEPREF]
                    [--callGeno CALLGENO] [--unif UNIF]
                    [--lrtSnp LRTSNP] [--lrtBia LRTBIA]
                    [--lrtTria LRTTRIA] [--minMaf MINMAF]
                    [--minQ MINQ]
                    [--minNonMajorCount MINNONMAJORCOUNT]
                    [--minNonMajorProportion MINNONMAJORPROPORTION]
                    [--minGlobalDepth MINGLOBALDEPTH]
                    [--maxGlobalDepth MAXGLOBALDEPTH]
                    [--minSampleDepth MINSAMPLEDEPTH]
                    [--maxSampleDepth MAXSAMPLEDEPTH]
                    [--minSamples MINSAMPLES] [--nGrids NGRIDS]
                    [--tol TOL] [--phredscale PHREDSCALE]
                    [--verbose VERBOSE] [--debug DEBUG]
                    [--printSites PRINTSITES] [-h]



optional arguments:
  --fin FIN             input file gzipped mpileup
  --fpars FPARS         pars (pANC and geno priors) file (default:
                        "NULL")
  --fout FOUT           output file gzipped text (default: "NULL")
  --fglikes FGLIKES     genotype likelihoods file gzipped text
                        (default: "NULL")
  --nSamples NSAMPLES   number of samples (type: Int64)
  --ploidy PLOIDY       ploidies to be tested, e.g. 2-5 or 1,3-5
                        (default: "1-8")
  --keepRef KEEPREF     keep reference as one possible allele (type:
                        Int64, default: 0)
  --callGeno CALLGENO   call genotypes (highest likelihood/posterior
                        probability) (type: Int64, default: 0)
  --unif UNIF           use a uniform prior for genotype probabilities
                        (type: Int64, default: 0)
  --lrtSnp LRTSNP       chisquare for SNP calling (type: Float64,
                        default: -Inf)
  --lrtBia LRTBIA       chisquare for biallelic calling (type:
                        Float64, default: -Inf)
  --lrtTria LRTTRIA     chisquare for triallelic (non) calling (type:
                        Float64, default: Inf)
  --minMaf MINMAF       minimum allele frequency for SNP calling
                        (type: Float64, default: -Inf)
  --minQ MINQ           minimum base quality in phredscore (type:
                        Int64, default: 13)
  --minNonMajorCount MINNONMAJORCOUNT
                        minimum non major base count (type: Int64,
                        default: 0)
  --minNonMajorProportion MINNONMAJORPROPORTION
                        minimum non major base proportion (type:
                        Float64, default: 0.0)
  --minGlobalDepth MINGLOBALDEPTH
                        minimum global depth (type: Int64, default: 1)
  --maxGlobalDepth MAXGLOBALDEPTH
                        maximum global depth (type: Int64, default:
                        10000)
  --minSampleDepth MINSAMPLEDEPTH
                        minimum sample depth (type: Int64, default: 0)
  --maxSampleDepth MAXSAMPLEDEPTH
                        maximum sample depth (type: Int64, default:
                        10000)
  --minSamples MINSAMPLES
                        minimum number of valid samples to retain site
                        (type: Int64, default: 1)
  --nGrids NGRIDS       grid density for grid-search estimation of
                        allele frequencies (type: Int64, default: 0)
  --tol TOL             tolerance for GSS estimation of allele
                        frequencies (type: Float64, default: 1.0e-5)
  --phredscale PHREDSCALE
                        phredscale (type: Int64, default: 33)
  --verbose VERBOSE     verbosity level (type: Int64, default: 0)
  --debug DEBUG         debug for 1 sample (type: Int64, default: 0)
  --printSites PRINTSITES
                        print on stdout every --printSites sites
                        (type: Int64, default: 10000)
  -h, --help            show this help message and exit
```






