
An alphabetical list of all available APIs is provided here.
Documentation can be assessed in `julia` using the help `?` command.

----------------------------------------------------------

- binomialExpansion(ploidy::Array{Int64,1}, freq::Float64)

Return the prior probabilities under HWE for a given ploidy and allele frequency.

```julia-repl
julia> binomialExpansion([2,3], 0.5) for diploid and triploid with allele frequency 0f 0.50
2-element Vector{Any}:
 Any[0.25, 0.5, 0.25]
 Any[0.125, 0.375, 0.375, 0.125]
```

---------------------------------------------------

- convertSyms(read::Reads, site::Site)

Convert symbols in the fifth element (sequence reads) of pileup (i.e. sequenced reads) to nucleotides from `Reads` and `Site` objects. Pileup format as defined here http://samtools.sourceforge.net/pileup.shtml

```julia-repl
rawReads=Reads(".G..G...","1533474323") # 10 reads and associated base qualities in Phred scores
mySite=Site("chrom12", 835132, 'A')
julia> convertSyms(rawReads, mySite)
("AGAAGAAA", Any[])
```

------------------------------------------------------------------------

- calcAlleleLike(read::Reads, allele::Array{Int64,1}; phredScale::Int64=33)

Calculate allele frequency likelihoods from a `Reads` object.

```julia-repl
julia> myReads=Reads("AGAAAGAAAA","1533474323")
julia> calcAlleleLike(myReads, [3], 1) # for allele G and haploid
-41.732245037944345
```

--------------------------------------------------------------------------

- calcFreqLike(read::Reads, allele::Array{Int64,1}, maf::Float64; phredScale::Int64=33)

Calculate the likelihood (in _ln_ format) for a given minor allele frequency, assuming haploid state.

```julia-repl
julia> myReads=Reads("AGAAAGAAAA","1533474323")
julia> calcFreqLike(myReads, [1,3], 0.10) # for major-minor alleles A and G with minor frequency 0.10
-5.545510518505054
```

----------------------------------------------------------------------------


- calcGenoLike(read::Reads, allele::Array{Int64,1}, ploidy::Int64; phredScale::Int64=33)

Calculate allele frequency likelihoods from a `Reads` object.

```julia-repl
julia> myReads=Reads("AGAAAGAAAA","1533474323")
julia> calcGenoLike(myReads, [1,3], 2) # for alleles A and G and diploid
3-element Vector{Float64}:
 -12.002917877610738
  -7.031999955835661
 -41.732245037944345
```

-----------------------------------------------------------------------------------

calcNonMajorCounts(read::Reads)

- Calculate the sum of non major alleles from a Reads object, useful to filter data based on the proportion (or count) of minor allele.

```julia-repl
julia> myReads=Reads("AGAAAGAAAA","1533474323")
julia> calcNonMajorCounts(myReads)
2
```

------------------------------------------------------------------------------------

- filterReads(read::Reads; phredScale::Int64=33, minBaseQuality::Int64=5)

Filter reads based on minimum base quality.

```julia-repl
julia> myReads=Reads("AGAAAGAAAA","1533474323")
julia> filterReads(myReads, minBaseQuality=20)
Reads("GG", "57")
```

-----------------------------------------------------------------------------------

- optimFreq(read::Reads, allele::Array{Int64,1}, tol::Float64)

Golden-search optimization for minor allele frequency given major and minor alleles.
Return likelihood value and the most likely minor allele frequency.

```julia-repl
julia> myReads=Reads("AGAAAGAAAA","1533474323")
julia> optimFreq(myReads, [1,3], 1e-4) # for major-minor alleles A and G with tolerance for low frequency of 1e-4
(-5.122068763452686, 0.19894055410912365)
```

-----------------------------------------------------------------------------------

- optimFreq\_GS(read::Reads, allele::Array{Int64,1}, nGrids::Int64)

Grid-search optimization for minor allele frequency given major and minor alleles.
Return likelihood value and the most likely minor allele frequency.

```julia-repl
julia> myReads=Reads("AGAAAGAAAA","1533474323")
julia> optimFreq_GS(myReads, [1,3], 100) # for major-minor alleles A and G with grid density of 1/100
(-5.122357762466471, 0.20202020202020202)
```

----------------------------------------------------------------------------

- snpTest(read::Reads, maxlike::Float64, allele::Array{Int64,1})

Return the (likelihood ratio test (LRT) statistic) for SNP calling.

```julia-repl
julia> myReads=Reads("AGAAAGAAAA","1533474323")
julia> freqsMLE=optimFreq(myReads, [1,3], 1e-4) # for alleles A and G
julia> snpTest(myReads, freqsMLE[1], [1, 3])
13.761698228316105
```

---------------------------------------------------------------------



