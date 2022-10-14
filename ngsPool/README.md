# ngsPool

Estimation of allele frequency (and more) from pooled sequencing data using `ngsJulia`.
`ngsPool` implements several estimators of allele frequencies from pooled NGS data.
It provides additional scripts for estimators of the site frequency spectrum (SFS) and for association tests.

To showcase its use, this tutorial will simulate some pooled sequencing data and demonstrate the various options and possible analyses implemented in `ngsPool`.
Before starting, please specify paths to both Julia language and ngsJulia (yours could be different):
```
JULIA=~/Software/julia-1.6.6/bin/julia
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
The file `test.txt` contains the true genotypes while `test.mpileup.gz` is a gzipped [mpileup](http://www.htslib.org/doc/samtools-mpileup.html) file containing information on sequencing data.
Specifically, columns in `text.txt` contain:
* contig identifier (set to the value of `--copy`)
* position
* reference allele (set to A)
* alternate allele (set to C)
* population allele frequency
* genotypes
* sample allele frequency

	
## Estimate allele frequencies with SNP calling and unknown sample size

Let's explore of `ngsPool` can be used to estimate allele frequencies.
We can retrieve a list of all available options by typing:
```
$JULIA $NGSJULIA/ngsPool/ngsPool.jl --help
```
The package requires a gzipped mpileup as input and the name of the output file in plain text format.
Several options for data filtering are available.
Let's understand its usage with several examples.

In many cases, the sample size is unknown. `ngsPool` provides a possibiity to obtain a maximum likelihood estimation (MLE) of the minor allele frequency under these circumstances.
Using the simulated data set, we can obtain per-site MLE of allele frequencies with:
```
$JULIA $NGSJULIA/ngsPool/ngsPool.jl --fin test.mpileup.gz --fout test.out.gz --lrtSnp 6.64
```
As an additional parameter, we specified a threshold for a Likelihood Ratio Test (LRT) for SNP calling.
In this example, the choice of 6.64 corresponds to a p-value of 0.01 (3.84 and 10.83 would correspond to p-values of 0.05 and 0.001, respectively).

The ouput file can be visualised with:
```
less -S test.out.gz
```
and for each called SNP provides the following information:
* chromosome
* position        
* reference allele
* nonreference allele
* major allele (inferred)
* minor allele (inferred) 
* lrtSNP (LRT statistic for SNP calling)
* lrtBia  (LRT statistic for bialleic site calling)
* lrtTria ((LRT statistic for trialleic site calling) 
* maf (estimated minor allele frequency)

The remaining columns are disabled using these options.

## Estimate allele frequency without SNP calling and known sample size

With known sample size, `ngsPool` calculates per-site sample allele frequency likelihoods which can be used to provide estimators of allele frequency or for further downstream analyses.
To this aim, the option `--nChroms` should be set equal to the product between ploidy and number of analysed samples.
Additionally if the option `--fsaf` is set, a gzipped text file with sample allele frequency likelihoods is returned.

Following our previous example, we can estimate allele frequencies (and their likelihoods) with:
```
$JULIA $NGSJULIA/ngsPool/ngsPool.jl --fin test.mpileup.gz --fout test.out.gz --nChroms 20 --fsaf test.saf.gz
```
Note that we are not doing SNP calling as ``--lrtSnp`` is not set.

The output file reports values for all sites that passed filtering:
```
less -S test.out.gz
```
and contains two additiona columns:
* saf_MLE (MLE of allele frequency from sample allele frequency likelihoods)
* saf_E (expected value of allele frequency from sample allele frequency likelihoods and uniform prior probability)

Additionally, a new file is generated:
```
less -S test.saf.gz
```
reporting the sample allele frequency log-likelihoods at each site (scaled to the ML).


## Site frequency spectrum (SFS)

The file containing the sample allele frequency log-likelihoods can be exploited for further downstream analyses.
For instance, `ngsPool` provides a script to estimate the SFS with three different methods.
This can be achieved with:
```
Rscript $NGSJULIA/ngsPool/poolSFS.R test.saf.gz > sfs.txt
```

The output file is accessible with
```
cat sfs.txt
```
and reports the estimated SFS based on:
* count: counting over MLE of per-site allele frequencies
* fit_count: fitting an exponential curve with counts of MLE of per-site allele frequencies
* fit_saf: fitting an exponential curve with per-site sample allele frequency likelihoods

# Association test

`ngsPool` provides a script to calculate association tests from sample allele frequency likelihoods.
Let's assume we have one target SNP and two groups, cases and controls, and we wish to test for a significant difference in allele frequencies.
We simulate different allele frequencis in two groups from low-depth pooled NGS data with
```
# cases
Rscript $NGSJULIA/simulMpileup_qq.R --out /dev/null --copy 2x200 --sites 1 --depth 1 --qq 0.1 --pool | gzip > test.cases.mpileup.gz

# controls
Rscript $NGSJULIA/simulMpileup_qq.R --out /dev/null --copy 2x200 --sites 1 --depth 1 --qq 0.05 --pool | gzip > test.controls.mpileup.gz
```

We calculate sample allele frequency likelihoods with:
```
# cases
$JULIA $NGSJULIA/ngsPool/ngsPool.jl --fin test.cases.mpileup.gz --fout /dev/null --nChroms 300 --fsaf test.cases.saf.gz 2> /dev/null

# controls
$JULIA $NGSJULIA/ngsPool/ngsPool.jl --fin test.controls.mpileup.gz --fout /dev/null --nChroms 300 --fsaf test.controls.saf.gz 2> /dev/null
```

These files are then used to test for association:
```
Rscript $NGSJULIA/ngsPool/poolAssoc.R test.cases.saf.gz test.controls.saf.gz > assoc.txt
```
and the resulting file is accessible with
```
cat assoc.txt
```
and shows LRT statistic and p-value (in log scale).
With multiple SNPs, each test result will be shown on different lines.






