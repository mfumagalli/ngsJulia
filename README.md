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
* [ngsPoly](https://github.com/mfumagalli/ngsJulia/tree/master/ngsPoly) infers the ploidy of samples from genotype likelihoods.
* [ngsPools](https://github.com/mfumagalli/ngsJulia/tree/master/ngsPool) estimates allele frequencies (and more) from pooled-sequencing data.





















