
# calculate geno likes and other stats from mpileup

include("../templates.jl")
include("../generics.jl")
include("functions.jl")

alleles = ['A','C','G','T']

# parsing parameters
parsed_args = parse_commandline_polyLite()
ploidy = parsePloidy(parsed_args["ploidy"])

## TEMPORARY, FIX THIS
if ploidy != [1,2,3,4,5,6,7,8]
	println("Sorry, the option --ploidy has to be fixed to 1-8 by the user at the moment.")
	exit(1)
end

MAXPLOIDY = 8
if maximum(ploidy) > MAXPLOIDY
	println("Error: tested ploidy must be less than", MAXPLOIDY+1)
	exit(1)
end

# open geno likes output file if set
if parsed_args["fglikes"] != "NULL"
        fglikes = GZip.open(parsed_args["fglikes"], "w")
end

# read mpileup
GZip.open(parsed_args["fin"]) do file

	for line in eachline(file)

		l = (split(line, "\t"))

		if parsed_args["nSamples"]!=((length(l)-3)/3)
			println("Error!: number of samples does not match input file.")
			println("Input:",(length(l)-3)/3)
			exit(1)
		end

		mySite = Site(l[1], parse(Int64, l[2]), uppercase(Char(l[3][1])))

		# pooled reads for first level filtering (global depth) and estimation of minor/major alleles and allele frequencies
		myReads = Reads("","")
		for n = 1:parsed_args["nSamples"]
			subReads = Reads(chomp(l[(n-1)*3+5]), chomp(l[(n-1)*3+6]))
			myReads.base = string(myReads.base, subReads.base)
			myReads.baseQuality = string(myReads.baseQuality, subReads.baseQuality)

		end

		# convert to bases
		(bases, indexDelN) = convertSyms(myReads, mySite)
		myReads = Reads(bases, myReads.baseQuality)

		error = 0
	  	# check if conversion was successful, if not print to stdout
		if ( !( length(myReads.base)==length(myReads.baseQuality)))
			println("Error! conversion:", myReads, mySite, "\n", length(myReads.base),"/",length(myReads.baseQuality))
			error = 1
		end

		# filter by quality
		if (error == 0)
			newReads = filterReads(myReads; phredScale=parsed_args["phredscale"], minBaseQuality=parsed_args["minQ"])
			myReads = newReads
			newReads = 0
		end
		globalDepth = length(myReads.base)

		# check if filtering was successful, if not print to stdout
		if ( !( length(myReads.base)==length(myReads.baseQuality)))
			if (error == 0)
				println("Error! filtering:", myReads, mySite, "\n", length(myReads.base),"/",length(myReads.baseQuality))
			end
			error = 1
		end

		# counts of non-major bases		
		nonMajorCount = calcNonMajorCounts(myReads)
		nonMajorProp = nonMajorCount/length(myReads.base)

		# filter the site based on global depth
		if (error == 0 && globalDepth>=parsed_args["minGlobalDepth"] && globalDepth<=parsed_args["maxGlobalDepth"] && (nonMajorCount>=parsed_args["minNonMajorCount"] || nonMajorCount==0) && (nonMajorProp>=parsed_args["minNonMajorProportion"] || nonMajorProp==0))

			# calculate genotype likelihoods assuming haploid for pooled samples
			haploid = calcGenoLogLike1(myReads, mySite)

			# order alleles
			if parsed_args["keepRef"] <= 0
				# calculate major and minor
				(major, minor, minor2, minor3) = sortperm(haploid, rev=true)
			else
				# reference is ancestral, set as major
				tmp_haploid = haploid
				tmp_haploid[(1:length(alleles))[mySite.reference .== alleles]]=0
				(major, minor, minor2, minor3) = sortperm(tmp_haploid, rev=true)
				tmp_haploid = 0
			end

			for n = 1:parsed_args["nSamples"]

				# retrieve bases for this particular sample
				myReads = Reads(chomp(l[(n-1)*3+5]), chomp(l[(n-1)*3+6]))

				# convert to bases
                		(bases, indexDelN) = convertSyms(myReads, mySite)
		       	        myReads = Reads(bases, myReads.baseQuality)

				suberror = 0
				# check if conversion was successful, if not print to stdout
				if ( !( length(myReads.base)==length(myReads.baseQuality)))
					println("Error! conversion:", myReads, mySite, "; sample:", n,"\n", length(myReads.base),"/",length(myReads.baseQuality))
					suberror = 1
				end

				# filter by quality
				if (suberror == 0)
					newReads = filterReads(myReads; phredScale=parsed_args["phredscale"], minBaseQuality=parsed_args["minQ"])
					myReads = newReads
					newReads = 0
				end
				sampleDepth = length(myReads.base)

				# check if filtering was successful, if not print to stdout
				if ( !( length(myReads.base)==length(myReads.baseQuality)))
					if (suberror == 0)
						println("Error! filtering:", myReads, mySite, "; sample:", n, "\n", length(myReads.base),"/",length(myReads.baseQuality))
					end
					suberror = 1
				end

				# haploid genotype likelihoods
				if 1 in ploidy
					haploid = calcGenoLogLike1_MajorMinor(myReads, mySite, major, minor)
				end

				# diploid genotype likelihoods
				if 2 in ploidy
					diploid = calcGenoLogLike2_MajorMinor(myReads, mySite, major, minor)
				end

				# triploid genotype likelihoods
				if 3 in ploidy
					triploid = calcGenoLogLike3_MajorMinor(myReads, mySite, major, minor)
				end

				# tetraploid genotype likelihoods
				if 4 in ploidy
					tetraploid = calcGenoLogLike4_MajorMinor(myReads, mySite, major, minor)
				end

				# pentaploid genotype likelihoods
				if 5 in ploidy
					pentaploid = calcGenoLogLike5_MajorMinor(myReads, mySite, major, minor)
				end

				# hexaploid genotype likelihoods
				if 6 in ploidy
					hexaploid = calcGenoLogLike6_MajorMinor(myReads, mySite, major, minor)
				end

				# heptaploid genotype likelihoods
                          	if 7 in ploidy
            	          		heptaploid = calcGenoLogLike7_MajorMinor(myReads, mySite, major, minor)
                            	end
					
				# octaploid genotype likelihoods
 				if 8 in ploidy
                                	octaploid = calcGenoLogLike8_MajorMinor(myReads, mySite, major, minor)
                             	end

				# write genotype likelihoods for each sample, if set
				if parsed_args["fglikes"]!="NULL"
					write(fglikes, join( (mySite.chrom, mySite.position, n, mySite.reference, sampleDepth, alleles[major], alleles[minor], join(haploid, "\t"), join(diploid, "\t"), join(triploid, "\t"), join(tetraploid, "\t"), join(pentaploid, "\t"), join(hexaploid, "\t"), join(heptaploid, "\t"), join(octaploid, "\t") ), "\t"), "\n")
				end

			end # for samples

		end # if not filtered for global depth

	end # for lines

end

if parsed_args["fglikes"] != "NULL"
        close(fglikes)
end


