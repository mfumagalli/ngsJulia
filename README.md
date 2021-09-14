# ngsJulia

Templates and functions in Julia language to process next-generation sequencing (NGS) data for population genetic analysis.
`ngsJulias` receives NGS data files as input and provides routines and functions to parse files, perform data filtering and implement custom population genetic analyses.
Two implementations for analysing pooled sequencing data and polyploid genomes are further presented.

## Installation

You can clone the repository with:
```
git clone https://github.com/mfumagalli/ngsJulia.git
```

If you wish to make sure you are using the most updated version you can do that with:
```
cd ngsJulia
git pull
```


## Dependencies

`ngsJulia` has been tested with Julia Version 1.6.1 (2021-04-23) available [here](https://julialang.org/downloads/).
It requires the packages `GZip` and `ArgParse` which can be obtained with:
```
julia> using Pkg
julia> Pkg.add("GZip")
julia> Pkg.add("ArgParse")
```

## Applications

We provide two novel applications of `ngsJulia` for low-coverage short-read sequencing data.
* [ngsPloidy](https://github.com/mfumagalli/ngsJulia/tree/master/ngsPloidy) infers the ploidy of samples from genotype likelihoods.
* [ngsPool](https://github.com/mfumagalli/ngsJulia/tree/master/ngsPool) estimates allele frequencies (and more) from pooled-sequencing data.

### Custom applications

`ngsJulia` has templates and functions that can be used to create custom analysis. 
As an illustration, assume we have sequencing data of a diallelic site for a __triploid__ organism and we wish to do genotype calling. 
Here how we can do it in `ngsJulia`.

Open a terminal in Julia (e.g., typing `~/Software/julia-1.6.1/bin/julia`) and load templates and functions in `ngsJulia`:
```
julia> include("templates.jl");
julia> include("functions.jl");
```

Let's assume we have the following sequencing data stored in these variables:
```
julia> myReads=Reads("AGAAAGAAAA","1533474323") # 10 reads and associated base qualities in Phred scores
julia> mySite=Site("chrom12", 835132, 'A') # chromosome, position and reference allele
```
These variable can be easily created by reading mpileup files, for instance using the following routine for this example:
```
using GZip

GZip.open("input.mpileup.gz") do file
	for line in eachline(file)
		l = (split(line, "\t"))
		global mySite = Site(l[1], parse(Int64, l[2]), uppercase(Char(l[3][1])))
		global myReads = Reads(chomp(l[5]), chomp(l[6]))
	end
end
```

We can visualise the nucleotide likelihoods:
```
julia> nucleo_likes = calcGenoLogLike1(myReads, mySite)
```
which in turn can be used to estimate major and minor alleles:
```
julia> (major, minor, minor2, minor3) = sortperm(nucleo_likes, rev=true);
julia> println("Major allele is ", alleles[major], " and minor allele is ", alleles[minor])
```

From these variables, it's easy to visualise the genotype likelihoods of a triploid for said alleles 
```
julia> geno_likes = calcGenoLogLike3_MajorMinor(myReads, mySite, major, minor)
```
where the genotypes in output are ordered as "(major,major,major), (major, major, minor), (major, minor, minor), (minor, minor, minor)", as that the most likely genotype is
```
julia> findmax(geno_likes)[2] # (major, major, minor), AAG
```

If we wish to set up a custom algorithm for genotype calling, then we can for instance calculate the difference in log likelihoods between the most likely and second most likely genotype as a weight of evidence
```
julia> diff(geno_likes[sortperm(geno_likes, rev=true)[[2,1]]])
```

In general, `ngsJulia` provides templates and functions useful for:
* data filtering based on quality and depth
* SNP and genotype calling
* nucleotide and genotype likelihoods for arbitrary ploidy
* allele frequency estimation

More specific usages can be found by investigating the code within `ngsPool` and `ngsPloidy`.

