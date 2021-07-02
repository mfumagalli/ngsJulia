
Things to add/do:

* in the (long) future this prototype will be replaced by a C++ software reading BAM files using bamtools

* fix the issue that so far only ploidy 1-6 is acceptable and valid

* add more possible ploidies, ideally in an automatic way, so up to an arbitrary number

* parallel implementation to make it faster

* test of issues on filtering based on global and sample depth, this has not been properly tested yet so bugs might be there

* do a lot of simulations to test its performance under various scenarios of ploidy and sequencing experiments; compute confusion matrices

* application on some real data

* add deviation from HWE option (impose inbreeding coefficients in writePars.R)

* joint estimation of ploidy and k (SFS shape) and F (inbreeding), EM algorithm?

* HMM to scan the chromosome and infer blocks of ploidy? Or simple sliding windows scan using whole-chrom likelihoods as prior
