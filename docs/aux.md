
Auxiliary `R` scripts provide additional functionalities for running simulations or CLIs, as listed here.

----------------------------------------------------------------------------------

- simulMpileup.R

Simulate mpileup file with different ploidy level for biallelic sites.

```bash
$ Rscript simulMpileup.R --help

out	o	2	character	output files for real data (and log if verbose), (mpileup is in stdout)
copy	c	1	character	ploidy per sample, e.g. 2x3,4 is 2,2,2,4
sites	s	2	integer	number of sites [default 1,000]
depth	d	2	double	mean haploid depth per sample [default 20.0]
lendepth	l	2	integer	mean length of sites with increasing/decreasing depth [default 0, disabled]
errdepth	e	2	double	error rate in mean depth [default 0.05]
qual	q	2	integer	mean base quality in phred score [default 20]
pvar	r	2	double	probability that site is variable in the population [1.0]
ksfs	k	2	double	coeff. for shape of SFS default [1.0]
panc	a	2	double	probability that ancestor state is correct [1.0]
ne	n	2	integer	effective population size [default 10,000]
pool	p	0	logical	enable pool data
help	h	0	logical	print help message
verbose	v	0	logical	verbose creates log file
offset	f	0	integer	offset value for genomic position
seed	u	2	integer	random seed for simulations reproducibility [default 180218]
```

```bash
$ Rscript simulMpileup.R --copy 3x2 --sites 3 --depth 5 # 2 triploid samples for 3 SNPs and average depth of 5X
copy_3x2	1	A	15	...............	543660320413024	17	.................	25234321684322154
copy_3x2	2	A	11	C...C.C.CC.	20441234234	19	C.C.CCCCCC.CC..CCCC	3262424423354344545
copy_3x2	3	A	13	........G....	4622432455514	11	.......CC..	42234625112
```


----------------------------------------------------------------------------------------

- writePars.R

In case of limited sample size, we need a create a file containing genotype probabilities.
We also need to provide the probability of the major allele being ancestral, as this information will be used in case of limited sample size.
These probability files can be generated using the following R script (which requires `getopt` package):

```bash
$ Rscript ngsPloidy/writePars.R --help
ksfs	k	2	double	coeff. for shape of SFS default [1.0]
ne	n	2	integer	effective population size [default 10,000]
snpcall	s	0	logical	flag id snps were called
panc	p	2	double	prob of ancestral being correct, if 0.5 is folded [default], if<0 compute
help	h	0	logical	print help message
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


```bash
$ Rscript ngsPloidy/writePars.R --ne 1000 --panc 0.5 # with effective population size of 1000 and folded spectrum
0.5	0.5
0.8665236	0.1334764
0.7508632	0.2313209	0.01781594
0.6506407	0.3006675	0.0463138	0.002378007
0.5637955	0.3473806	0.08026401	0.008242398	0.0003174078
0.4885422	0.3762669	0.1159178	0.01785558	0.001375207	0.00004236644
0.4233333	0.391253	0.1506682	0.03094457	0.003574947	0.0002202691	0.000005654918
0.3668283	0.395535	0.1827806	0.04692484	0.007228144	0.0006680394	0.00003430084	0.0000007547979
0.3178654	0.3917033	0.2111783	0.06505838	0.01252672	0.001543658	0.00011889	0.000005232402	0.0000001007477
```

Further examples are:

* with known allelic polarisation:
```bash
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 1
```
* with uncertain assignment of ancestral alleles (probability of being correct of 0.90)
```bash
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 0.90
```
* with automatic calculation of probability of misassignment of ancestral allele
```bash
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p -1
```
* with folded spectrum (unknown allelic polarisation)
```bash
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 0.5
```
* or with SNPs only
```bash
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 1 -s
```

-------------------------------------------------------------------------------------------------




