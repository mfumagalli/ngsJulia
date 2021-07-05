# Estimate ploidy levels from gzipped mpileup files

# run 'julia ngsPloidy --help' for documentation
include("../structure.jl")
include("../functions.jl")
include("args.jl")

# alleles = ['A','C','G','T']

# parsing parameters
parsed_args = parse_commandline_poly()
ploidy = parsePloidy(parsed_args["ploidy"])

# at the moment only 1-8 ploidies are accepted
if ploidy != [1,2,3,4,5,6,7,8]
	println("Sorry, the option --ploidy has to be fixed to 1-8 by the user at the moment.")
	exit(1)
end

MAXPLOIDY = 8
if maximum(ploidy) > MAXPLOIDY
	println("Error: tested ploidy must be less than", MAXPLOIDY+1)
	exit(1)
end

# check that if you given genotype priors, then you assume that reference in mpileup is ancestral
if parsed_args["fpars"] != "NULL" && parsed_args["keepRef"] <= 0
	println("Error: if you give genotype priors with --pars, then --keepRef must be 1.")
	exit(1)
end
# likewise if you don't give genotype priors, then you don't care about anc/derived but sample size must be high
if parsed_args["fpars"] == "NULL" && parsed_args["nSamples"] < 10
	if parsed_args["debug"] != 2
		println("Error: if you don't give genotype priors with --pars, then the recommended minimum sample size is 10.")
		exit(1)
	end
end

# if you use a uniform prior, then you need all samples to have data, otherwise missing data will results in equally probable genotypes
if parsed_args["unif"] > 0 && parsed_args["minSamples"] < parsed_args["nSamples"]
	println("Error: if you use a uniform prior for genotype probabilities then you must require that all samples have data, therefore you should set --minSamples equal to --nSamples.")
	exit(1)
end
# note that in case of equally probable genotypes and genotype calling, the homozygous ancestral/reference/most-common-allele will be chosen

# note on lrt values:
# TH = 7.879 is pv = 0.005
# TH = 6.635 is pv = 0.01
# TH = 3.841 is pv = 0.05
# TH = -1 means disabled (any negative number)

# read priors on genotypes and probability that major is ancestral
pANC = 1.0 # if no geno priors are given, then no need to switch anc/der
if parsed_args["fpars"]!="NULL"
	(pANC, priors) = readPriors(parsed_args["fpars"], ploidy)
end

# print on screen
if parsed_args["verbose"] > 0
	println("Parsed args:")
	for (arg,val) in parsed_args
		println("  $arg  =>  $val")
	end
	println(ploidy)
	if parsed_args["fpars"] != "NULL"
		println(pANC, "\n", priors)
	end
end

# open output file if set
if parsed_args["fout"] != "NULL"
	fout = GZip.open(parsed_args["fout"], "w")
	if parsed_args["keepRef"] <= 0
		write(fout, join( ("##chrom", "pos", "ref", "depth", "major", "minor", "lrtSnp", "lrtBia", "lrtTria", "maf"), "\t"), "\n")
	else
		write(fout, join( ("##chrom", "pos", "ref", "depth", "ref", "alt", "lrtSnp", "lrtBia", "lrtTria", "aaf"), "\t"), "\n")
	end
end

# open geno likes output file if set
if parsed_args["fglikes"] != "NULL"
        fglikes = GZip.open(parsed_args["fglikes"], "w")
end

# maximum not-valid samples
maxNotPassed = parsed_args["nSamples"] - parsed_args["minSamples"]

# initialise the nr of informative (used) sites
infSites = 0
# initialise matrix of likelihoods: nr_ploidy^nsamples possible events
# matrix samples on rows, ploidies on columns
polyLikes = zeros(parsed_args["nSamples"], length(ploidy))

# read mpileup
GZip.open(parsed_args["fin"]) do file

	infSites = zeros(parsed_args["nSamples"])
	allSites = 0
	polyLikes = zeros(parsed_args["nSamples"], length(ploidy))

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

			allSites += 1

			# calculate genotype likelihoods assuming haploid for pooled samples
			haploid = calcGenoLogLike1(myReads, mySite)

			# order alleles
			if parsed_args["keepRef"] <= 0
				# calculate major and minor
				(major, minor, minor2, minor3) = sortperm(haploid, rev=true)
			else
				# reference is ancestral, set as major
				tmp_haploid = haploid
			  	tmp_haploid[findmax(mySite.reference .== alleles)[2]] = 0 # set it max
				(major, minor, minor2, minor3) = sortperm(tmp_haploid, rev=true)
				tmp_haploid = 0
			end

			# is this biallelic?
			if parsed_args["thBia"] > -Inf
				biaLike = calcGenoLogLike1_Bia(myReads, mySite, major, minor, phredScale=parsed_args["phredscale"])
				# if you fix the ref/anc, then the "minor" may have greater likelihood, so:
				if parsed_args["keepRef"] <= 0
					lrtBia = 2*(biaLike-haploid[major])
				else
					lrtBia = maximum(2*(biaLike-haploid[major]), 2*(haploid[major]-biaLike))
				end
			else
				lrtBia = Inf
			end

			# is this triallelic?
			if parsed_args["thTria"] < Inf
				triaLike = calcGenoLogLike1_Tria(myReads, mySite, major, minor, minor2, phredScale=parsed_args["phredscale"])
				lrtTria = 2*(triaLike-biaLike)
			else
				lrtTria = -Inf
			end

			# is it polymorphic?
			if parsed_args["nGrids"] > 0
				freqsMLE = optimFreq_MajorMinor(myReads, major, minor, parsed_args["nGrids"]+1)
			else
				freqsMLE = optimFreq_MajorMinor_GSS(myReads, major, minor, parsed_args["tol"])
			end
			lrtSnp = snpPval_MajorMinor(myReads, freqsMLE[1], major, minor)

			# SNP calling
			if (lrtSnp>parsed_args["thSnp"] && lrtBia>parsed_args["thBia"] && lrtTria<parsed_args["thTria"] && freqsMLE[2]>=parsed_args["minMaf"] && (1-freqsMLE[2])>=parsed_args["minMaf"]  )

				# initialise iterator to check that all (or the required number of) samples pass filtering
				samplesPassed = 0
				# initiate the per-site likelihoods, rows are samples, columns are ploidies
				polySite = zeros(parsed_args["nSamples"], length(ploidy))
				# is it a SNP within each sample?
				snpSite = zeros(parsed_args["nSamples"])

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

					# skip site if errors in one sample
					if suberror == 1
						n = parse_args["nSamples"] + 1 # terminate the loop
					end

					# increment nr of valid samples
					if (suberror == 0 && sampleDepth >= parsed_args["minSampleDepth"] && sampleDepth <= parsed_args["maxSampleDepth"])
						samplesPassed += 1
					end

					# to save computational time, skip if already enough samples do not pass filtering
					if (suberror) == 0 && (n - samplesPassed <= maxNotPassed)

						snpSite[n] = 1

						# haploid genotype likelihoods
						if 1 in ploidy
							haploid = calcGenoLogLike1_MajorMinor(myReads, mySite, major, minor)
							if maximum(pANC) < 1
								haploid_rev = calcGenoLogLike1_MajorMinor(myReads, mySite, minor, major)
							else
								haploid_rev = zeros(length(haploid))
							end
						end

						# diploid genotype likelihoods
						if 2 in ploidy
							diploid = calcGenoLogLike2_MajorMinor(myReads, mySite, major, minor)
							if maximum(pANC) < 1
								diploid_rev = calcGenoLogLike2_MajorMinor(myReads, mySite, minor, major)
							else
								diploid_rev = zeros(length(diploid))
							end
						end

						# triploid genotype likelihoods
						if 3 in ploidy
							triploid = calcGenoLogLike3_MajorMinor(myReads, mySite, major, minor)
							if maximum(pANC) < 1
								triploid_rev = calcGenoLogLike3_MajorMinor(myReads, mySite, minor, major)
							else
								triploid_rev = zeros(length(triploid))
							end
						end

						# tetraploid genotype likelihoods
						if 4 in ploidy
							tetraploid = calcGenoLogLike4_MajorMinor(myReads, mySite, major, minor)
							if maximum(pANC) < 1
								tetraploid_rev = calcGenoLogLike4_MajorMinor(myReads, mySite, minor, major)
							else
								tetraploid_rev = zeros(length(tetraploid))
							end
						end

						# pentaploid genotype likelihoods
						if 5 in ploidy
							pentaploid = calcGenoLogLike5_MajorMinor(myReads, mySite, major, minor)
							if maximum(pANC) < 1
								pentaploid_rev = calcGenoLogLike5_MajorMinor(myReads, mySite, minor, major)
							else
								pentaploid_rev = zeros(length(pentaploid))
							end
						end

						# hexaploid genotype likelihoods
						if 6 in ploidy
							hexaploid = calcGenoLogLike6_MajorMinor(myReads, mySite, major, minor)
							if maximum(pANC) < 1
								hexaploid_rev = calcGenoLogLike6_MajorMinor(myReads, mySite, minor, major)
							else
								hexaploid_rev = zeros(length(hexaploid))
							end
						end

						# heptaploid genotype likelihoods
                                                if 7 in ploidy
                                                        heptaploid = calcGenoLogLike7_MajorMinor(myReads, mySite, major, minor)
                                                        if maximum(pANC) < 1
                                                                heptaploid_rev = calcGenoLogLike7_MajorMinor(myReads, mySite, minor, major)
                                                        else
                                                                heptaploid_rev = zeros(length(heptaploid))
                                                        end
                                                end
					
						# octaploid genotype likelihoods
                                                if 8 in ploidy
                                                        octaploid = calcGenoLogLike8_MajorMinor(myReads, mySite, major, minor)
                                                        if maximum(pANC) < 1
                                                                octaploid_rev = calcGenoLogLike8_MajorMinor(myReads, mySite, minor, major)
                                                        else
                                                                octaploid_rev = zeros(length(octaploid))
                                                        end
                                                end

						# allele frequency (either derived or alternate or minor)
						freq =  freqsMLE[2]

						for p in 1:length(ploidy)

							# which ploidy
							genoLikes = genoLikes_rev = Array{Float64}
							if ploidy[p] == 1
								genoLikes = haploid
								genoLikes_rev = haploid_rev
								genoPriors = (1-freq,freq)
							elseif ploidy[p] == 2
								genoLikes = diploid
								genoLikes_rev = diploid_rev
								genoPriors = ( (1-freq)^2, 2*(1-freq)*freq, freq^2 )
							elseif ploidy[p] == 3
								genoLikes = triploid
								genoLikes_rev = triploid_rev
								genoPriors = ( (1-freq)^3, 3*(1-freq)^2*freq, 3*(1-freq)*freq^2, freq^3 )
							elseif ploidy[p] == 4
								genoLikes = tetraploid
								genoLikes_rev = tetraploid_rev
								genoPriors = ( (1-freq)^4, 4*(1-freq)^3*freq, 6*(1-freq)^2*freq^2, 4*(1-freq)*freq^3, freq^4 )
							elseif ploidy[p] == 5
								genoLikes = pentaploid
								genoLikes_rev = pentaploid_rev
								genoPriors = ( (1-freq)^5, 5*(1-freq)^4*freq, 10*(1-freq)^3*freq^2, 10*(1-freq)^2*freq^3, 5*(1-freq)*freq^4, freq^5 )
							elseif ploidy[p] == 6
								genoLikes = hexaploid
								genoLikes_rev = hexaploid_rev
								genoPriors = ( (1-freq)^6, 6*(1-freq)^5*freq, 15*(1-freq)^4*freq^2, 20*(1-freq)^3*freq^3, 15*(1-freq)^2*freq^4, 6*(1-freq)*freq^5, freq^6 )
							elseif ploidy[p] == 7
                                                                genoLikes = heptaploid
                                                                genoLikes_rev = heptaploid_rev
                                                                genoPriors = ( (1-freq)^7, 7*(1-freq)^6*freq, 21*(1-freq)^5*freq^2, 35*(1-freq)^4*freq^3, 35*(1-freq)^3*freq^4, 21*(1-freq)^2*freq^5, 7*(1-freq)*freq^6, freq^7 )

							elseif ploidy[p] == 8
                                                        	genoLikes = octaploid
                                                      		genoLikes_rev = octaploid_rev
                                                     		genoPriors = ( (1-freq)^8, 8*(1-freq)^7*freq, 28*(1-freq)^6*freq^2, 56*(1-freq)^5*freq^3, 70*(1-freq)^4*freq^4, 56*(1-freq)^3*freq^5, 28*(1-freq)^2*freq^6, 8*(1-freq)*freq^7, freq^8 )


							end # match ploidy

							# change prior if uniform is set
							if parsed_args["unif"]==1
								genoPriors = ones(Float64, ploidy[p]+1) / (ploidy[p]+1)
							end

							# calculate all probabilities
							tempProbs = []
							for z = 1:length(genoLikes)
								if parsed_args["fpars"] != "NULL"
									tempProbs = [tempProbs; exp(genoLikes[z]) * priors[p,z] * pANC[1]; exp(genoLikes_rev[z]) * priors[p,z] * pANC[2] ]
	              else
									tempProbs = [tempProbs; exp(genoLikes[z]) * genoPriors[z] ]
								end
							end

							# assign genotypes or integrate across all of them
							if parsed_args["callGeno"] > 0
								polySite[n,p] += maximum(tempProbs)
							else
								polySite[n,p] += sum(tempProbs)
							end

						end # loop ploidy

						# print debug
						if parsed_args["debug"]==1
							println(mySite, myReads, ";", log(polySite))
						end
	
						# write genotype likelihoods for each sample, if set
						if parsed_args["fglikes"]!="NULL"
							write(fglikes, join( (mySite.chrom, mySite.position, n, mySite.reference, sampleDepth, alleles[major], alleles[minor], join(haploid, "\t"), join(diploid, "\t"), join(triploid, "\t"), join(tetraploid, "\t"), join(pentaploid, "\t"), join(hexaploid, "\t"), join(heptaploid, "\t"), join(octaploid, "\t") ), "\t"), "\n")
						end

					else # if not enough samples

						n = parsed_args["nSamples"] + 1 # terminate the loop

					end # if enough or not samples

				end # for samples

				# test if enough samples passed the filter and if so update the polyLikes matrix
				if samplesPassed >= parsed_args["minSamples"]

					infSites += snpSite
					for i in 1:size(polySite,1)
						for j in 1:size(polySite,2)
							if polySite[i,j]!=0
								# sum across site
								polyLikes[i,j] += log(polySite[i,j])
							end
						end
					end

					if parsed_args["verbose"]>1
						println(l, "\n", log(polySite))
					end

					if parsed_args["fout"]!="NULL"
						if parsed_args["nSamples"] > 1
							write(fout, join( (mySite.chrom, mySite.position, mySite.reference, globalDepth, alleles[major], alleles[minor], lrtSnp, lrtBia, lrtTria, freqsMLE[2]), "\t"), "\n")
						else
							write(fout, join( (mySite.chrom, mySite.position, mySite.reference, globalDepth, alleles[major], alleles[minor], lrtSnp, lrtBia, lrtTria, freqsMLE[2], log.(polySite[1]), log.(polySite[2]), log.(polySite[3]), log.(polySite[4]), log.(polySite[5]), log.(polySite[6]), log.(polySite[7]), log.(polySite[8])  ), "\t"), "\n")
						end

					end

          				if parsed_args["verbose"]>0 && (allSites % parsed_args["printSites"]) == 0
						println("allSites:", allSites, "; infSites:", infSites, "; last:", mySite.position)
						println(polyLikes)
						for p in 1:length(ploidy)
							println("LIKE all ploidy ", ploidy[p], ":", sum(polyLikes[:,p]) )
						end
					end

				end # if completed

			end # if SNP

		end # if not filtered for global depth

	end # for lines

end

if parsed_args["fout"] != "NULL"
	close(fout)
end

if parsed_args["fglikes"] != "NULL"
        close(fglikes)
end

if parsed_args["debug"]<2
	println("infSites:", infSites)
end

# fix NaN to -Inf
for i in 1:size(polyLikes,1)
	for j in 1:size(polyLikes,2)
		if isnan(polyLikes[i,j])
			polyLikes[i,j]=-Inf
		end
	end
end

println("LIKELIHOODS:\n", polyLikes)

# do not calculate AIC/BIC anymore?
#(aic, bic) = calcModelStat(polyLikes, ploidy, infSites)

#println("AIC:\n", aic)
#println("BIC:\n", bic)

# maximum likelihood
maxLike = 0.0
samPloidy = zeros(parsed_args["nSamples"])
for n in 1:parsed_args["nSamples"]
	samPloidy[n] = sortperm(vec(polyLikes[n,:]), rev=true)[1]
	maxLike = maxLike + polyLikes[n,sortperm(vec(polyLikes[n,:]), rev=true)[1]]
end

println("Max like ploidy:", samPloidy, "\nwith log-like:", maxLike)

# all ploidy - maximum likelihood
for p in 1:length(ploidy)
	println("(LIKE-max like) all ploidy ", ploidy[p], ":", sum(polyLikes[:,p])-maxLike )
end



