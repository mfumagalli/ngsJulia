
`ngsJulia` has templates and functions that can be used to create custom analysis.
As an illustration, assume we have sequencing data of a diallelic site for a __triploid__ organism and we wish to do genotype calling.
Here how we can do it in `ngsJulia`.

Open a Julia'shell and load `ngsJulia`:
```julia-repl
julia> include("ngsJulia.jl");
```

Let's assume we have the following sequencing data stored in these variables:
```julia-repl
julia> myReads=Reads("AGAAAGAAAA","1533474323"); # 10 reads and associated base qualities in Phred scores
julia> mySite=Site("chrom12", 835132, 'A'); # chromosome, position and reference allele
```
These variable can be easily created by reading mpileup files, for instance using the following routine for this example:
```julia-repl
>julia 	GZip.open("input.mpileup.gz") do file
        	for line in eachline(file)
                	l = (split(line, "\t"))
                	global mySite = Site(l[1], parse(Int64, l[2]), uppercase(Char(l[3][1])))
                	global myReads = Reads(chomp(l[5]), chomp(l[6]))
        	end
	end
```

We can visualise the nucleotide likelihoods:
```julia-repl
julia> nucleoLikes = [calcGenoLike(myReads, [i], 1) for i=1:4];
```
which in turn can be used to estimate major and minor alleles:
```julia-repl
julia> (major, minor, minor2, minor3) = sortperm(nucleoLikes, rev=true);
julia> println("Major allele is ", ALLELES[major], " and minor allele is ", ALLELES[minor]);
```

From these variables, it's easy to visualise the genotype likelihoods of a triploid for said alleles
```julia-repl
julia> genoLikes = calcGenoLike(myReads, [major, minor], 3);
```
where the genotypes in output are ordered as "(major,major,major), (major, major, minor), (major, minor, minor), (minor, minor, minor)", as that the most likely genotype is
```julia-repl
julia> findmax(genoLikes)[2] # (major, major, minor), AAG
```

If we wish to set up a custom algorithm for genotype calling, then we can for instance calculate the difference in log likelihoods between the most likely and second most likely genotype as a weight of evidence
```julia-repl
julia> diff(genoLikes[sortperm(genoLikes, rev=true)[[2,1]]])
```

In general, `ngsJulia` provides templates and functions useful for:
* data filtering based on quality and depth
* SNP and genotype calling
* nucleotide and genotype likelihoods for arbitrary ploidy
* allele frequency estimation




