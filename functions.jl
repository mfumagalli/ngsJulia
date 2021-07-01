
# Generic Functions

# convert ascii character to phred score (this function is deprecated, and not used)
#function ascii2phred(ascii::Char)

#	asciiSymbols = """!"#\$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~""";

#	return search(asciiSymbols, ascii)

#end

# calculate genotype likelihoods (in ln format) in case of haploids
function calcGenoLogLike1(read::Reads, site::Site, phredScale::Int64=33)
	#the only function not requiring major and minor, as major and minor
	 #allele is defined based on the haploid genotype likelihood as output of it

	alleles = ['A', 'C', 'G', 'T']
	# initialise likelihoods
	likes = zeros(length(alleles))

	# cycle across all possible genotypes
	for j1 = 1:length(alleles)

		# cycle across all reads
		for i = 1:length(read.base)

			# get base probability from quality score
			bP = 10^((phredScale - Int64(read.baseQuality[i])  )/10) # mapping quality scaled in the phred scale to obtain the error probability of wrong alignment
			# mapping quality of a read alignment: (phredScale - Int64(read.baseQuality[i])
			# p: p of that alignment is wrong

			sublike = 0.0
			if alleles[j1]==read.base[i]
				sublike += 1-(bP)
			else
				sublike += (bP/3)
			end
			likes[j1] += log(sublike)
		end
	end
	return likes #genotype likelihood vector
end


# calculate genotype likelihoods (in ln format) in case of haploids for
 # EITHER major OR minor (for SNP calling, biallelic test)
function calcGenoLogLike1_Bia(read::Reads, site::Site, major::Int64, minor::Int64; phredScale::Int64=33)

	alleles = ['A', 'C', 'G', 'T']
	like = 0.0

	# cycle across all reads
	for i = 1:length(read.base)
		# get base probability from quality score
		bP = 10^((phredScale - Int64(read.baseQuality[i]))/10)
		sublike = 0.0;
		if read.base[i]==alleles[major] || read.base[i]==alleles[minor]
			sublike += 1-(bP)
		else
			sublike += (bP/3)
		end
		like += log(sublike)
	end
	return like
end

# calculate genotype likelihoods (in ln format) in case of haploids for
 # EITHER major OR minor OR second minor (for SNP calling, triallelic test)
function calcGenoLogLike1_Tria(read::Reads, site::Site, major::Int64, minor::Int64, minor2::Int64; phredScale::Int64=33)

	alleles = ['A', 'C', 'G', 'T']
	like = 0.0

	# cycle across all reads
	for i = 1:length(read.base)
		# get base probability from quality score
		bP = 10^((phredScale - Int64(read.baseQuality[i]))/10)
		sublike = 0.0
		if read.base[i]==alleles[major] || read.base[i]==alleles[minor] || read.base[i]==alleles[minor2]
			sublike += 1-(bP)
		else
			sublike += (bP/3)
		end
		like += log(sublike)
	end
	return like
end


# calculate genotype likelihoods (in ln format) in case of haploids for
 # Major and Minor
function calcGenoLogLike1_MajorMinor(read::Reads, site::Site, major::Int64, minor::Int64, phredScale::Int64=33)

	alleles = ['A', 'C', 'G', 'T']
	likes = zeros(2)
	iter = 0
	ploidy = 1

	# cycle across all possible genotypes
	for (j1) = ((major), (minor))
		iter +=1 ;
	        # cycle across all reads
		for i = 1:length(read.base)
			# get base probability from quality score
			bP = 10^((phredScale - Int64(read.baseQuality[i]))/10)
			sublike = 0.0
			if alleles[j1]==read.base[i]
				sublike += (1-bP)/ploidy
			else
				sublike += (bP/3)/ploidy
			end
			likes[iter] += log(sublike)
		end
	end
	return likes
end;


# calculate genotype likelihoods (in ln format) in case of diploids for Major and Minor
function calcGenoLogLike2_MajorMinor(read::Reads, site::Site, major::Int64, minor::Int64, phredScale::Int64=33)

	likes = zeros(3)
	alleles = ['A', 'C', 'G', 'T']
	iter = 0
	ploidy = 2

	# cycle across only considered genotypes
	for (j1,j2) = ((major,major),(major,minor),(minor,minor))

		iter += 1
		# cycle across all reads
		for i = 1:length(read.base)

			# calculate base probability
	        	bP = 10^((phredScale - Int64(read.baseQuality[i]))/10)

			sublike = 0.0
	        if alleles[j1]==read.base[i]
				sublike += (1-bP)/ploidy
			else
				sublike += (bP/3)/ploidy
			end
			if alleles[j2]==read.base[i]
				sublike += (1-bP)/ploidy
			else
				sublike += (bP/3)/ploidy
			end
			likes[iter] += log(sublike)
		end
	end
	return likes
end;

# calculate genotype likelihoods (in ln format) in case of triploids for Major and Minor
function calcGenoLogLike3_MajorMinor(read::Reads, site::Site, major::Int64, minor::Int64, phredScale::Int64=33)

	likes = zeros(4);
	alleles = ['A', 'C', 'G', 'T'];
	iter = 0;
	ploidy = 3

	# cycle across only considered genotypes
	for (j1,j2,j3) = ((major,major,major),(major,major,minor),(major,minor,minor),(minor,minor,minor))

		iter += 1;
		# cycle across all reads
		for i = 1:length(read.base)

			bP = 10^((phredScale - Int64(read.baseQuality[i]))/10);

	                sublike = 0.0;
			if alleles[j1]==read.base[i]
				sublike += (1-bP)/ploidy;
			else
				sublike += (bP/3)/ploidy;
			end
			if alleles[j2]==read.base[i]
				sublike += (1-bP)/ploidy;
			else
				sublike += (bP/3)/ploidy;
			end
			if alleles[j3]==read.base[i]
				sublike += (1-bP)/ploidy;
		    	else
				sublike += (bP/3)/ploidy;
			end
	 		likes[iter] += log(sublike);
		end
	end

	return likes
end;

# calculate genotype likelihoods (in ln format) in case of tetraploids for Major and Minor
function calcGenoLogLike4_MajorMinor(read::Reads, site::Site, major::Int64, minor::Int64, phredScale::Int64=33)

	likes = zeros(5);
	alleles = ['A','C','G','T'];
	iter = 0;
	ploidy = 4;

	# cycle across only considered genotypes
	for (j1,j2,j3,j4) = ((major,major,major,major),(major,major,major,minor),(major,major,minor,minor),(major,minor,minor,minor),(minor,minor,minor,minor))

		iter += 1;
		for i = 1:length(read.base)

			bP = 10^((phredScale - Int64(read.baseQuality[i]))/10);

			sublike = 0.0;
			if alleles[j1]==read.base[i]
				sublike += (1-bP)/ploidy;
			else
				sublike += (bP/3)/ploidy;
			end
			if alleles[j2]==read.base[i]
				sublike += (1-bP)/ploidy;
			else
				sublike += (bP/3)/ploidy;
			end
			if alleles[j3]==read.base[i]
				sublike += (1-bP)/ploidy;
			else
				sublike += (bP/3)/ploidy;
			end
			if alleles[j4]==read.base[i]
				sublike += (1-bP)/ploidy;
			else
				sublike += (bP/3)/ploidy;
			end
			likes[iter] += log(sublike);
		end
	end
	return likes;
end;

# # calculate genotype likelihoods (in ln format) in case of pentaploids for Major and Minor
function calcGenoLogLike5_MajorMinor(read::Reads, site::Site, major::Int64, minor::Int64, phredScale::Int64=33)

	likes = zeros(6);
	alleles = ['A','C','G','T'];
	iter = 0;
	ploidy = 5

	# cycle across only considered genotypes
        for (j1,j2,j3,j4,j5) = ((major,major,major,major,major),(major,major,major,major,minor),(major,major,major,minor,minor),(major,major,minor,minor,minor),(major,minor,minor,minor,minor),(minor,minor,minor,minor,minor))

		iter += 1;
		for i = 1:length(read.base)

			bP = 10^((phredScale - Int64(read.baseQuality[i]))/10);

			sublike = 0.0;
			if alleles[j1]==read.base[i]
				sublike += (1-bP)/ploidy;
			else
				sublike += (bP/3)/ploidy;
			end
			if alleles[j2]==read.base[i]
				sublike += (1-bP)/ploidy;
			else
				sublike += (bP/3)/ploidy;
			end
			if alleles[j3]==read.base[i]
				sublike += (1-bP)/ploidy;
			else
				sublike += (bP/3)/ploidy;
			end
			if alleles[j4]==read.base[i]
				sublike += (1-bP)/ploidy;
			else
				sublike += (bP/3)/ploidy;
			end
			if alleles[j5]==read.base[i]
				sublike += (1-bP)/ploidy;
			else
				sublike += (bP/3)/ploidy;
			end
			likes[iter] += log(sublike);
		end
	end
	return likes;
end;

# # calculate genotype likelihoods (in ln format) in case of sessaploids for Major and Minor
function calcGenoLogLike6_MajorMinor(read::Reads, site::Site, major::Int64, minor::Int64, phredScale::Int64=33)

	likes = zeros(7);
	alleles = ['A','C','G','T'];
	iter = 0;
	ploidy = 6

	# cycle across only considered genotypes
	for (j1,j2,j3,j4,j5,j6) = ((major,major,major,major,major,major),(major,major,major,major,major,minor),(major,major,major,major,minor,minor),(major,major,major,minor,minor,minor),(major,major,minor,minor,minor,minor),(major,minor,minor,minor,minor,minor),(minor,minor,minor,minor,minor,minor))

		iter += 1;
		for i = 1:length(read.base)

			bP = 10^((phredScale - Int64(read.baseQuality[i]))/10);
			sublike = 0.0;

			if alleles[j1]==read.base[i]
				sublike += (1-bP)/ploidy
			else
				sublike += (bP/3)/ploidy
			end
			if alleles[j2]==read.base[i]
				sublike += (1-bP)/ploidy
			else
				sublike += (bP/3)/ploidy
			end
			if alleles[j3]==read.base[i]
				sublike += (1-bP)/ploidy
			else
				sublike += (bP/3)/ploidy
			end
			if alleles[j4]==read.base[i]
				sublike += (1-bP)/ploidy
			else
				sublike += (bP/3)/ploidy
			end
			if alleles[j5]==read.base[i]
				sublike += (1-bP)/ploidy
			else
				sublike += (bP/3)/ploidy
			end
			if alleles[j6]==read.base[i]
				sublike += (1-bP)/ploidy
			else
				sublike += (bP/3)/ploidy
			end
			likes[iter] += log(sublike)
		end
	end
	return likes
end;

# # calculate genotype likelihoods (in ln format) in case of settaploids for Major and Minor
function calcGenoLogLike7_MajorMinor(read::Reads, site::Site, major::Int64, minor::Int64, phredScale::Int64=33)

        likes = zeros(8)
        alleles = ['A','C','G','T']
        iter = 0
	ploidy = 7

	# cycle across only considered genotypes
        for (j1,j2,j3,j4,j5,j6,j7) = ((major,major,major,major,major,major,major),(major,major,major,major,major,major,minor),(major,major,major,major,major,minor,minor),(major,major,major,major,minor,minor,minor),(major,major,major,minor,minor,minor,minor),(major,major,minor,minor,minor,minor,minor),(major,minor,minor,minor,minor,minor,minor), (minor,minor,minor,minor,minor,minor,minor) )

                iter += 1;
                for i = 1:length(read.base)

                        bP = 10^((phredScale - Int64(read.baseQuality[i]))/10);
                        sublike = 0.0;

                        if alleles[j1]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else                                                                                                                             			sublike += (bP/3)/ploidy
                        end
                        if alleles[j2]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
                        if alleles[j3]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
                        if alleles[j4]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
                        if alleles[j5]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
                        if alleles[j6]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
			if alleles[j7]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
                        likes[iter] += log(sublike)
                end
        end
	return likes
end

# # calculate genotype likelihoods (in ln format) in case of octaploids for Major and Minor
function calcGenoLogLike8_MajorMinor(read::Reads, site::Site, major::Int64, minor::Int64, phredScale::Int64=33)

        likes = zeros(9)
        alleles = ['A','C','G','T']
        iter = 0
        ploidy = 8

        # cycle across only considered genotypes
        for (j1,j2,j3,j4,j5,j6,j7,j8) = ((major,major,major,major,major,major,major,major),(major,major,major,major,major,major,major,minor),(major,major,major,major,major,major,minor,minor),(major,major,major,major,major,minor,minor,minor),(major,major,major,major,minor,minor,minor,minor),(major,major,major,minor,minor,minor,minor,minor),(major,major,minor,minor,minor,minor,minor,minor), (major,minor,minor,minor,minor,minor,minor,minor), (minor,minor,minor,minor,minor,minor,minor,minor) )

                iter += 1;
                for i = 1:length(read.base)

                        bP = 10^((phredScale - Int64(read.baseQuality[i]))/10);
                        sublike = 0.0;

                        if alleles[j1]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
							sublike += (bP/3)/ploidy
                        end
                        if alleles[j2]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
                        if alleles[j3]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
                        if alleles[j4]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
                        if alleles[j5]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
						if alleles[j6]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
                        if alleles[j7]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
						if alleles[j8]==read.base[i]
                                sublike += (1-bP)/ploidy
                        else
                                sublike += (bP/3)/ploidy
                        end
                        likes[iter] += log(sublike)
                end
        end
        return likes
end


# convert symbols in fifth element(sequence reads) in pileup to nucleotides
# pileup format as defined here: http://samtools.sourceforge.net/pileup.shtml
function convertSyms(read::Reads, site::Site)

	bases = ""
	indexDelN = Array{Int64, 1} #position on original pileup read

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
			if parse(Int, read.base[i+2])<10 # indel longer than 9
				println("generic_line485  ", typeof(read.base[i+1]))
				lenIndel = parse(Int, string(read.base[i+1], read.base[i+2])) + 1 #Matteo: digit, so +1 #?
			end
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

# filter based on proportion or count of minor allele (major substracted)
function calcNonMajorCounts(read::Reads)

	alleles = ['A','C','G','T']
	counts = Int.(zeros(4))

	if length(read.base)>0
		for i = 1:length(read.base)
			counts[alleles.==split(read.base,"")[i][1]] += 1
		end
	end

	# return sum(counts)-counts[sortperm(counts, rev=true)[1]]
	return sum(counts)-sort!(counts)[1]
end


# calculate likelihood (in ln format) for a given minor allele frequency
 # (in case of haploids) from major and minor
function calcFreqLogLike1_MajorMinor(read::Reads, major::Int64, minor::Int64, maf::Float64, phredScale::Int64=33)
	#in calcFreqLogLike1_MajorMinor() maf input set to 0.0

	like = 0
	alleles = ['A','C','G','T']
	freqs = [1-maf, maf] #major, minor allele frequency

	for i = 1:length(read.base)
		bP::Float64 = 10^(  (phredScale - Int64(read.baseQuality[i])  )/10)
		sublike = 0.0
		iter = 0
		for j = (major,minor)
			iter += 1
			if alleles[j]==read.base[i]
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
function optimFreq_MajorMinor(read::Reads, major::Int64, minor::Int64, nGrids::Int64)

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
		like=calcFreqLogLike1_MajorMinor(read, major, minor, maf)
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
function optimFreq_MajorMinor_GSS(read::Reads, major::Int64, minor::Int64, tol::Float64) #tol: tolerance

	maxlike = -Inf
	MLEmaf = Float64

	a = 0.0 # minimum
	b = 1.00 # maximum

	# golden ratio
	gr = (sqrt(5)+1)/2

	c = b -  (b - a) / gr
	d = a + (b - a) / gr

	while abs(c-d) > tol

		fc = calcFreqLogLike1_MajorMinor(read, major, minor, c) #c and d as maf
		fd = calcFreqLogLike1_MajorMinor(read, major, minor, d)

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

	maxlike = calcFreqLogLike1_MajorMinor(read, major, minor, MLEmaf)

	return (maxlike, MLEmaf)

end

# LRT(likelihood ratio test statistic) for snp calling
function snpPval_MajorMinor(read::Reads, maxlike::Float64, major::Int64, minor::Int64)

	maf = 0.0
	like_H0 = maximum( [calcFreqLogLike1_MajorMinor(read, major, minor, maf);
	calcFreqLogLike1_MajorMinor(read, minor, major, maf)] )
	# global and sample major could be differ

	lrtSNP =  2 * (maxlike - like_H0 ) #-2(ln(p_H0))-ln(parameter space)) #LRT statistic
	#freqMLE=maxlike=ln(parameter space)

	return lrtSNP
end
