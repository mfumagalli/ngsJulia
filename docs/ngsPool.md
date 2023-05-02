# ngsPool

`ngsPool` implements a method to estimate allele frequencies and perform various analyses from pooled sequencing data using `ngsJulia`.
The model is described in the related [paper](https://f1000research.com/articles/11-126).
In addition to provide several estimators of allele frequencies, `ngsPool` includes scripts to estimate the site frequency spectrum and for association tests.

To showcase its use, this tutorial will simulate some pooled sequencing data and demonstrate the various options and possible analyses implemented in `ngsPool`.
Throughout these examples, we assume that we defined an environment variable `NGSJULIA` that points to the installation path.


## Simulate pooled NGS data 

We can simulate NGS data from a pooled sequencing experiment using a script provided in `ngsJulia`.
We can explore its options:
```bash
Rscript $NGSJULIA/simulMpileup.R --help
```
which are also accessible on this [page](https://ngsjulia.readthedocs.io/en/latest/aux/).

Let's assume that we wish to simulate 10 diploid genomes, 1000 base pairs each with an average sequencing depth of 20 and base quality of 20 in Phred score. Samples come from a population of 10,000 effective size under constant-size evolution.
We can do that by running:
```bash
Rscript $NGSJULIA/simulMpileup.R --out test.txt --copy 2x10 --sites 1000 --depth 20 --qual 20 --ksfs 1 --ne 10000 --pool | gzip > test.mpileup.gz
```
and explore the files generated with:
```bash
ls test.*
```
The file `test.txt` contains the true genotypes while `test.mpileup.gz` is a gzipped [mpileup](http://www.htslib.org/doc/samtools-mpileup.html) file containing information on sequencing data.
Specifically, columns in `text.txt` contain:  
 contig identifier (set to the value of `--copy`)  
-position  
-reference allele (set to A)  
-alternate allele (set to C)  
-population allele frequency  
-genotypes  
-sample allele frequency

	
## Estimate allele frequencies with SNP calling and unknown sample size

Let's explore of `ngsPool` can be used to estimate allele frequencies.
We can retrieve a list of all available options by typing:
```bash
julia $NGSJULIA/ngsPool/ngsPool.jl --help
```
The package requires a gzipped mpileup as input and the name of the output file in plain text format.
Several options for data filtering are available.
Let's understand its usage with several examples.

In many cases, the sample size is unknown. `ngsPool` provides a possibiity to obtain a maximum likelihood estimation (MLE) of the minor allele frequency under these circumstances.
Using the simulated data set, we can obtain per-site MLE of allele frequencies with:
```bash
julia $NGSJULIA/ngsPool/ngsPool.jl --fin test.mpileup.gz --fout test.out.gz --lrtSnp 6.64
```
As an additional parameter, we specified a threshold for a Likelihood Ratio Test (LRT) for SNP calling.
In this example, the choice of 6.64 corresponds to a p-value of 0.01 (3.84 and 10.83 would correspond to p-values of 0.05 and 0.001, respectively).

The ouput file can be visualised with:
```bash
less -S test.out.gz
```
and for each called SNP provides the following information:  
-chromosome
-position      
-reference allele  
-nonreference allele  
-major allele (inferred)  
-minor allele (inferred)  
-lrtSNP (LRT statistic for SNP calling)  
-lrtBia  (LRT statistic for bialleic site calling)  
-lrtTria ((LRT statistic for trialleic site calling)   
-maf (estimated minor allele frequency)

The remaining columns are disabled using these options.

## Estimate allele frequency without SNP calling and known sample size

With known sample size, `ngsPool` calculates per-site sample allele frequency likelihoods which can be used to provide estimators of allele frequency or for further downstream analyses.
To this aim, the option `--nChroms` should be set equal to the product between ploidy and number of analysed samples.
Additionally if the option `--fsaf` is set, a gzipped text file with sample allele frequency likelihoods is returned.

Following our previous example, we can estimate allele frequencies (and their likelihoods) with:
```bash
julia $NGSJULIA/ngsPool/ngsPool.jl --fin test.mpileup.gz --fout test.out.gz --nChroms 20 --fsaf test.saf.gz
```
Note that we are not doing SNP calling as ``--lrtSnp`` is not set.

The output file reports values for all sites that passed filtering:
```bash
less -S test.out.gz
```
and contains two additiona columns:
* saf\_MLE (MLE of allele frequency from sample allele frequency likelihoods)
* saf\_E (expected value of allele frequency from sample allele frequency likelihoods and uniform prior probability)

Additionally, a new file is generated:
```bash
less -S test.saf.gz
```
reporting the sample allele frequency log-likelihoods at each site (scaled to the ML).

## Site frequency spectrum

The file containing the sample allele frequency log-likelihoods can be exploited for further downstream analyses.
For instance, `ngsPool` provides a script to estimate the SFS with three different methods.
This can be achieved with:
```bash
Rscript $NGSJULIA/ngsPool/poolSFS.R test.saf.gz > sfs.txt
```

The output file is accessible with
```bash
cat sfs.txt
```
and reports the estimated SFS based on:  
-count: counting over MLE of per-site allele frequencies  
-fit\_count: fitting an exponential curve with counts of MLE of per-site allele frequencies  
-fit\_saf: fitting an exponential curve with per-site sample allele frequency likelihoods

## Association test

`ngsPool` provides a script to calculate association tests from sample allele frequency likelihoods.
Let's assume we have one target SNP and two groups, cases and controls, and we wish to test for a significant difference in allele frequencies.
We simulate different allele frequencis in two groups from low-depth pooled NGS data with
```bash
# cases
Rscript $NGSJULIA/simulMpileup_qq.R --out /dev/null --copy 2x200 --sites 1 --depth 1 --qq 0.1 --pool | gzip > test.cases.mpileup.gz

# controls
Rscript $NGSJULIA/simulMpileup_qq.R --out /dev/null --copy 2x200 --sites 1 --depth 1 --qq 0.05 --pool | gzip > test.controls.mpileup.gz
```

We calculate sample allele frequency likelihoods with:
```bash
# cases
julia $NGSJULIA/ngsPool/ngsPool.jl --fin test.cases.mpileup.gz --fout /dev/null --nChroms 300 --fsaf test.cases.saf.gz 2> /dev/null

# controls
julia $NGSJULIA/ngsPool/ngsPool.jl --fin test.controls.mpileup.gz --fout /dev/null --nChroms 300 --fsaf test.controls.saf.gz 2> /dev/null
```

These files are then used to test for association:
```bash
Rscript $NGSJULIA/ngsPool/poolAssoc.R test.cases.saf.gz test.controls.saf.gz > assoc.txt
```
and the resulting file is accessible with
```bash
cat assoc.txt
```
and shows LRT statistic and p-value (in log scale).
With multiple SNPs, each test result will be shown on different lines.

## Further options

All options available in `ngsPool` can be accessed with:
```bash
julia ngsPool/ngsPool.jl --help

usage: ngsPool.jl --fin FIN --fout FOUT [--fsaf FSAF]
                  [--nChroms NCHROMS] [--lrtSnp LRTSNP]
                  [--lrtBia LRTBIA] [--lrtTria LRTTRIA] [--minQ MINQ]
                  [--minDepth MINDEPTH] [--maxDepth MAXDEPTH]
                  [--nGrids NGRIDS] [--tol TOL]
                  [--phredscale PHREDSCALE] [--verbose VERBOSE]
                  [--printSites PRINTSITES] [-h]

optional arguments:
  --fin FIN             input file gzipped mpileup
  --fout FOUT           output file gzipped text
  --fsaf FSAF           output gzipped saf file (default: "/dev/null")
  --nChroms NCHROMS     total number of chromosomes pooled (ploidy *
                        number of individuals) [>0 ensables saf
                        likelihoods] (type: Int64, default: 0)
  --lrtSnp LRTSNP       LRT for SNP calling (type: Float64, default:
                        -Inf)
  --lrtBia LRTBIA       LRT for biallelic calling (type: Float64,
                        default: -Inf)
  --lrtTria LRTTRIA     LRT for triallelic (non) calling (type:
                        Float64, default: Inf)
  --minQ MINQ           minimum base quality in phredscore (type:
                        Int64, default: 5)
  --minDepth MINDEPTH   minimum global depth (type: Int64, default: 1)
  --maxDepth MAXDEPTH   maximum global depth (type: Int64, default:
                        100000)
  --nGrids NGRIDS       grid density for grid-search estimation of
                        allele frequencies (type: Int64, default: 0)
  --tol TOL             tolerance for GSS estimation of allele
                        frequencies (type: Float64, default: 1.0e-5)
  --phredscale PHREDSCALE
                        phredscale (type: Int64, default: 33)
  --verbose VERBOSE     verbosity level (type: Int64, default: 1)
  --printSites PRINTSITES
                        print on stdout every --printSites sites
                        (type: Int64, default: 10000)
  -h, --help            show this help message and exit

```


