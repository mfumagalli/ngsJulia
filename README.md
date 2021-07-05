# ngsJulia

Templates and functions in Julia language to process next-generation sequencing data for population genetic analysis.

## Installation

```
git clone https://github.com/mfumagalli/ngsJulia.git
```

## Dependencies

`ngsJulia` has been tested with Julia Version 1.6.1 (2021-04-23) available [here](https://julialang.org/downloads/).
It requires the packages `GZip` and `ArgParse` which can be obtained with:
```
julia> using Pkg
julia> Pkg.add("GZip")
julia> Pkg.add("GZip")
```

## Applications

We provide two novel applications of `ngsJulia` for low-coverage short-read sequencing data.
* [ngsPloidy](https://github.com/mfumagalli/ngsJulia/tree/master/ngsPloidy) infers the ploidy of samples from genotype likelihoods.
* [ngsPools](https://github.com/mfumagalli/ngsJulia/tree/master/ngsPool) estimates allele frequencies (and more) from pooled-sequencing data.

### Custom applications

`ngsJulia` has templates and functions that can be used to create custom analysis. For instance, imagine we have sequencing data of a diallelic site for a triploid organism and we wish to do genotype calling. Here how we can do it in `ngsJulia`.

```
# load templates and functions
julia> include("structure.jl")
julia> include("functions.jl")

# assume we the following sequencing data stored in these variables:
julia> myReads=Reads("AAGAGGAAAC","5342123560") # 10 reads and associated base qualities in Phred scores
julia> mySite=Site("chrom12", 34512, 'A') # chromosome, position and reference allele

# these variable can be easily created by reading mpileup files, for instance using the following routine for one sample:
#GZip.open(parsed_args["input.mpileup.gz"]) do file
#	for line in eachline(file)
#		l = (split(line, "\t"))
#		myReads = Site(l[1], parse(Int64, l[2]), uppercase(Char(l[3][1])))
#		mySite = Reads(chomp(l[5]), chomp(l[6]))
#	end
#end

# we can visualise the nucleotide likelihoods
julia> nucleo_likes = calcGenoLogLike1(myReads, mySite)

# which in turn can be used to estimate major and minor alleles
julia> (major, minor, minor2, minor3) = sortperm(nucleo_likes, rev=true)
julia> println("Major allele is ", alleles[major], " and minor allele is ", alleles[minor])

# from this, it's easy to visualise the genotype likelihoods of a triploid for said alleles 
julia> geno_likes = calcGenoLogLike3_MajorMinor(myReads, mySite, major, minor)
# where the genotypes in output are ordered as "(major,major,major), (major, major, minor), (major, minor, minor), (minor, minor, minor)"
# as such the most likely genotype is
julia> findmax(geno_likes)[2] # (major, major, minor), AAG

# if we wish to set up a custom algorithm for genotype calling, then we can for instance calculate the difference in log likelihoods between the most likely and second most likely genotype as a weight of evidence
julia> diff(geno_likes[sortperm(geno_likes, rev=true)[[2,1]]])
```

In general, `ngsJulia` provides templates and functions useful for:
* data filtering based on quality and depth
* SNP and genotype calling
* nucleotide and genotype likelihoods for arbitrary ploidy
* allele frequency estimation

More specific usages can be found by investigating the code within `ngsPool` and `ngsPloidy`.






















calcNonMajorCounts(myReads)







```



















