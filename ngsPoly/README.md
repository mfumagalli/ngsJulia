# ngsPoly
Inference of ploidy from short-read sequencing data

Collaborators, please check and edit my [WISHLIST](wishlist.md) and [RANDOM](random.md) thoughts.

## Model

Notation:
* P: probability
* D: data
* O: ploidy
* A: ancestral state
```
log(P(D, G, A|O_1=y_1, O_2=y_2, ..., O_n=y_n)) = sum_{samples} sum_{sites} P(D|G,A) P(G|O,A) P(A)
```
	
where:
* P(D|G,A): genotype likelihood (assuming two alleles at most)
* P(G|O,A): genotype probability (derived from the expected population allele frequency from the site frequency spectrum)
* P(A): probability of correct assignment of ancestral state

We assume we observed at most two alleles.

The program will output the most likely array of marginal ploidies, as well as the likelihood of all samples having the same ploidy.
In the latter assumption, we can propose a Bayesian formulation:

```
P(O_1=y, O_2=y, ..., O_n=y | D) = P(D|O) P(O) / P(D)
```

where P(D|O) is calculated as aforementioned.

## Installation

This has been tested with Julia >= 0.4.7 (Version 0.4.7-pre+3, Commit eb322c0\*) on a x86\_64-linux-gnu machine and it requires the following packages: GZip, ArgParse.

	git clone https://github.com/mfumagalli/ngsPoly.git
        git clone https://github.com/mfumagalli/ngsJulia.git

        cd ngsPoly
        ln -s ../ngsJulia/simulMpileup.R simulMpileup.R
        ln -s ../ngsJulia/generics.jl generics.jl
        ln -s ../ngsJulia/templates.jl templates.jl

## Example

First, you need a create a file containing your prior probabilities and genotypes and major allele being the ancestral.
This can be generated using the following R script.

	Rscript writePars.R -k 1 -n 10000 -p 1 > test.pars
	Rscript writePars.R -k 1 -n 10000 -p 0.90 > test.unk.pars
	Rscript writePars.R -k 1 -n 10000 -p -1 > test.auto.pars
	Rscript writePars.R -k 1 -n 10000 -p 0.5 > test.fold.pars
	Rscript writePars.R -k 0.9 -n 100000 -p 1 -s > test.snp.pars

where: 
* '-k' denotes the shape of the site frequency spectrum:
	- k=1 : constant population size
	- k>1 : population bottleneck
	- k<1 : population growth
* '-n' is the effective population size
* '-s' if a flag and specifies that SNPs have beeen called; actually this makes sense only if only one samples is analysed with called SNPs 
* '-p' specifies how to define the probability that the ancestral state is the major allele, if 0.5 this means you assume folded data, if -1 it will compute it using the site frequency spectrum
* '-h' prints a help message.

You can simulate NGS data in mpileup format by specifying the ploidy of each individual and other parameters of the sequencing experiment
and species.
Run `Rscript simulMpileup.R --help` for help.
After that, you can run `julia ngsPoly.jl --help` to see all options for estimating ploidy.

### Case A: 2 haploids, 2 diploids, 2 triploids, 2 tetraploids, 2 pentaploids

	Rscript simulMpileup.R --out test.A.txt --copy 1x2,2x2,3x2,4x2,5x2 --sites 5000 --depth 100 --qual 20 --ksfs 1 --ne 10000 | gzip > test.A.mpileup.gz

	less -S test.A.txt

	zcat test.A.mpileup.gz | less -S

	# known ancestral state (reference)
	julia ngsPoly.jl --fin test.A.mpileup.gz --fpars test.pars --fout test.A.out.gz --nSamples 10 --thSnp -1 --ploidy 1-5
	# automatic set of probability of major allele being the ancestral state
	julia ngsPoly.jl --fin test.A.mpileup.gz --fpars test.auto.pars --fout test.A.out.gz --nSamples 10 --thSnp -1 --ploidy 1-5 --keepRef 0

Please note that if `--fout` is given, the program will print some statistics for each site, including the estimate allele frequency.

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

People contacted for potential application on real data and/or sugegstions and advice (in no precise order):

* Simon O'Hanlon (Imperial), frog fungus
* Peter Fields (UoBasel), snails
* Melissa Wilson-Sayres (ASU), sex-linked
* Benjamin Schwessinger (ANU), plant-fungi
* Reuben Nowell (Silwood), bdelloid rotifers

People to contact:

* Yoshida et al. eLife 2013 (via Ben)
* Lien Bertier (UC Davis, via Ben)



