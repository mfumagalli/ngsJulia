# ngsJulia

Templates and functions in Julia language to process next-generation sequencing (NGS) data for population genetic analysis.
`ngsJulias` receives NGS data files as input and provides routines and functions to parse files, perform data filtering and implement custom population genetic analyses.
Two implementations for analysing pooled sequencing data and polyploid genomes are further presented.

Full documentation is accessible [here](https://ngsjulia.readthedocs.io).

You can clone the repository with `git clone https://github.com/mfumagalli/ngsJulia.git`

We also provide two novel applications of `ngsJulia` for low-coverage short-read sequencing data.
* [ngsPloidy](https://github.com/mfumagalli/ngsJulia/tree/master/ngsPloidy) infers the ploidy of samples from genotype likelihoods.
* [ngsPool](https://github.com/mfumagalli/ngsJulia/tree/master/ngsPool) estimates allele frequencies (and more) from pooled-sequencing data.

Archived code and scripts to replicate all results in the accompanying paper are available the `paper` folder.


