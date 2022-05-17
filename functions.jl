
function calcAlleleLike(read::Reads, allele::Array{Int64}; phredScale::Int64=33)

        like = 0.0

	for i = 1:length(read.base) # each based
        	bP = 10^((phredScale - Int64(read.baseQuality[i]))/10)
                sublike = 0.0
                if read.base[i] in ALLELES[allele]
                	sublike += (1-bP)/ploidy
            	else
                  	sublike += (bP/3)/ploidy
                end
                like += log(sublike)
        end
        return like
end

function calcGenoLike(read::Reads, allele::Array{Int64}, ploidy::Int64; phredScale::Int64=33)

        likes = zeros(ploidy+1)
        iter = 0
        genotypes = collect(with_replacement_combinations(allele, ploidy))

        for genotype in genotypes
                iter += 1
		for i = 1:length(read.base) # each based
                        bP = 10^((phredScale - Int64(read.baseQuality[i]))/10)
                        sublike = 0.0
                        for item in genotype # each base in genotype
                                if ALLELES[item]==read.base[i]
                                        sublike += (1-bP)/ploidy
                                else
                                        sublike += (bP/3)/ploidy
                                end
                        end
                        likes[iter] += log(sublike)
                  end
	end
	if (length(allele)==1) 
		likes = likes[1]
	end
        return likes
end

# convert symbols in fifth element(sequence reads) in pileup to nucleotides
# pileup format as defined here: http://samtools.sourceforge.net/pileup.shtml
function convertSyms(read::Reads, site::Site)

	bases = ""
	#indexDelN = Array{Int64, 1} #position on original pileup read
	indexDelN = []

	i=1
	while i <= length(read.base)

		if read.base[i] in ['.',','] # reference
			bases=string(bases,site.reference)
		elseif uppercase(read.base[i]) in ['A','C','G','T'] # alternate
			bases=string(bases, uppercase(read.base[i]))
		elseif read.base[i] in ['^']
		 # start read, index +1 since the following character is mapping quality
			i = i + 1
		# elseif read.base[i] in ['$'] # end read
		# 	i = i # do nothing but keep for clarity
		elseif read.base[i] in ['*', 'N', 'n', '>', '<']
			 # asterisk is deleted base, N or n is undefined, > and < are
			  # reference skips, these will be then filtered out later on
			bases = string(bases, "X")
			 # but keep track of this as there is associated base quality
			# indexDelN = [indexDelN; i]
			push!(indexDelN, i)
		elseif read.base[i] in ['-', '+']
		 # indel, skip to the the next non-indel base
			lenIndel = parse(Int, read.base[i+1])
		  	#if parse(Int, read.base[i+2])<10 # indel longer than 9
		  	#	println("generic_line485  ", typeof(read.base[i+1]))
		  	#	lenIndel = parse(Int, string(read.base[i+1], read.base[i+2])) + 1 #Matteo: digit, so +1 #?
		  	#end
			#if parse(Int, read.base[i+3])<10 # indel longer than 99
			#	lenIndel = parse(Int, string(read.base[i+1], read.base[i+2], read.base[i+3]))
			#end
			i = i + lenIndel + 1 #  no associated base quality
		end
		i += 1
	end

	return (bases, indexDelN)
end

# filter based on minimum base quality
function filterReads(read::Reads; phredScale::Int64=33, minBaseQuality::Int64=5)

	filt=Reads("","")

	if (length(read.base)>0)
		i=1;
		while i <= length(read.base)
		 # filtering for base quality (also check is it is not an X
		  # my notation for either deleted base or N)
       			if read.base!="X" &&  (Int64(read.baseQuality[i])-phredScale) >= minBaseQuality
				#The ASCII of the character following `^' minus 33 gives the mapping quality.
				#mapping quality >= minBaseQuality

				filt.base=string(filt.base, read.base[i])
				filt.baseQuality=string(filt.baseQuality, read.baseQuality[i])
			end
			i += 1
		end
	end
	return filt
end

"""
	calcNonMajorCounts(read::Reads)

This function receives a read object and returns the sum of non major alleles.
This is used to filter data based on the proportion (or count) of minor allele.
"""
function calcNonMajorCounts(read::Reads)

	counts = Int.(zeros(4))

	if length(read.base)>0
		for i = 1:length(read.base)
		     counts += (ALLELES.==split(read.base,"")[i])
		end
	end

	return sum(counts)-sort!(counts, rev=true)[1]
end


# calculate likelihood (in ln format) for a given minor allele frequency
 # (in case of haploids) from major and minor
function calcFreqLike(read::Reads, allele::Array{Int64}, maf::Float64; phredScale::Int64=33)
	#in calcFreqLogLike1_MajorMinor() maf input set to 0.0

	like = 0
	
	freqs = [1-maf, maf] #major, minor allele frequency

	for i = 1:length(read.base)
		bP::Float64 = 10^(  (phredScale - Int64(read.baseQuality[i])  )/10)
		sublike = 0.0
		iter = 0
		for j = allele
			iter += 1
			if ALLELES[j]==read.base[i]
				sublike += (1 - bP)*freqs[iter]
			else
				sublike += (bP/3)*freqs[iter]
			end
		end
		like += log(sublike)
	end

	return like
end


#suspected to be an unfinished function..(2020/5/27)
# grid-search optimization for allele frequencies, only for Major and Minor alleles
function optimFreq_GS(read::Reads, allele::Array{Int64}, nGrids::Int64)

	maxlike = -Inf
	# MLEmaf = Float64 #what's this for? collect? array or a number?
	#  how about a tuple?
	MLEmaf = Float64[] #is this necessary or not?

	# define the grid
	# fGrid = linspace(0,1.0,nGrids)
	fGrid = collect(range(0,1,length=nGrids))

	i = 1;
	# while fGrid[i] <= 1.00
	while i <= nGrids

		# maf = Float64(fGrid[i])
		maf = fGrid[i] #as already float64
		like=calcFreqLike(read, allele, maf)
		if like > maxlike #not filtering sites at all
			maxlike=like
			MLEmaf=maf #intending for the last grid num?
		end
		i += 1
	end
  # #the while (actually for loop) make no sense as it will always take maf=1
  # maf=1
  # MLEmaf=maf #intending for the last grid num?
  # maxlike = calcFreqLogLike1_MajorMinor(read, major, minor, maf)

	return (maxlike, MLEmaf)
end

# golden section search optimization for allele frequencies
# only for Major and Minor alleles
function optimFreq(read::Reads, allele::Array{Int64}, tol::Float64) #tol: tolerance

	maxlike = -Inf
	MLEmaf = Float64

	a = 0.0 # minimum
	b = 1.00 # maximum

	# golden ratio
	gr = (sqrt(5)+1)/2

	c = b -  (b - a) / gr
	d = a + (b - a) / gr

	while abs(c-d) > tol

		fc = calcFreqLike(read, allele, c) #c and d as maf
		fd = calcFreqLike(read, allele, d)

		if fc > fd # to find the maximum, otherwise use < (Matteo) ??
			b = d
		else #fc <= fd
			a = c
		end
		# we recompute both c and d here to avoid loss of precision which
		# may lead to incorrect results or infinite loop
		# as suggested on wiki page
		c = b - (b - a) / gr
		d = a + (b - a) / gr
	end

	MLEmaf = (b+a)/2

	# set to 0 if very low, below some threshold based on the tolerance
	if MLEmaf < tol*2
		MLEmaf = 0.0
	end

	maxlike = calcFreqLike(read, allele, MLEmaf)

	return (maxlike, MLEmaf)

end

# LRT(likelihood ratio test statistic) for snp calling
function snpTest(read::Reads, maxlike::Float64, allele::Array{Int64})

	maf = 0.0
	like_H0 = maximum( [calcFreqLike(read, [major, minor], maf)
	calcFreqLike(read, [minor, major], maf)] )
	# global and sample major could be differ

	lrtSNP =  2 * (maxlike - like_H0 ) #-2(ln(p_H0))-ln(parameter space)) #LRT statistic
	#freqMLE=maxlike=ln(parameter space)

	return lrtSNP
end



