
# SNP calling, allele frequency estimation, sample allele frequency likelihoods for pool sequencing data (from gzipped mpileup)

# run 'julia ngsPool --help' for documentation

include("templates.jl")
include("generics.jl")
include("functions.jl")

alleles = ['A','C','G','T']

# parsing parameters
parsed_args = parse_commandline_pool()
println("Parsed args:")
for (arg,val) in parsed_args
	println("  $arg  =>  $val")
end
if parsed_args["nChroms"] < 0
	parsed_args["nChroms"] = 0
end

# open output file for writing
f2 = GZip.open(parsed_args["fout"], "w")
# write header
write(f2, join( ("##chrom", "position", "reference", "nonreference", "major", "minor", "lrtSNP", "lrtBia", "lrtTria", "maf", "naf_max", "naf_E"), "\t"), "\n")

if parsed_args["fsaf"]!="/dev/null"
	f3 = GZip.open(parsed_args["fsaf"], "w")
end

# read mpileup
GZip.open(parsed_args["fin"]) do file

	allSites = infSites = 0

	for line in eachline(file)

		l = (split(line, "\t"))

		mySite = Site(l[1], parse(Int64, l[2]), Char(l[3][1]))
		myReads = Reads(l[5], chomp(l[6])) 
		# convert to bases
		(bases, indexDelN) = convertSyms(myReads, mySite);
		myReads = Reads(bases, myReads.baseQuality);
		# if there is a deleted base or a N, the base is marked as X and will be filtered out in the next stage

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
		depth = length(myReads.base)

		# check if filtering was successful, if not print to stdout
		if ( !( length(myReads.base)==length(myReads.baseQuality)))
			if (error == 0) 
				println("Error! filtering:", myReads, mySite, "\n", length(myReads.base),"/",length(myReads.baseQuality))
			end
			error = 1
		end

		# filter by depth
		if (error == 0 && depth >= parsed_args["minDepth"] && depth <= parsed_args["maxDepth"])

			# valid sites
			allSites += 1

			# calculate genotype likelihoods for all haploid cases
			haploid=calcGenoLogLike1(myReads, mySite)
			# order alleles
			(major, minor, minor2, minor3) = sortperm(haploid, rev=true)
				
			# is this biallelic?
			biaLike = calcGenoLogLike1_Bia(myReads, mySite, major, minor, phredScale=parsed_args["phredscale"])
			lrtBia = 2*(biaLike-haploid[major])

			# is this triallelic?
			triaLike = calcGenoLogLike1_Tria(myReads, mySite, major, minor, minor2, phredScale=parsed_args["phredscale"])
			lrtTria = 2*(triaLike-biaLike)

			# is this a SNP?
			if parsed_args["nGrids"] > 0
				freqsMLE = optimFreq_MajorMinor(myReads, major, minor, parsed_args["nGrids"]+1)
			else
				freqsMLE = optimFreq_MajorMinor_GSS(myReads, major, minor, parsed_args["tol"])
			end
			lrtSNP = snpPval_MajorMinor(myReads, freqsMLE[1], major, minor)

			if (lrtSNP>parsed_args["thSnp"] && lrtBia>parsed_args["thBia"] && lrtTria<parsed_args["thTria"])

				# used sites
				infSites += 1
					
				# calculate saf (sample allele frequencies) likelihoods for the non-ref allele
				freqMax = -1.0
				freqE = -1.0
				safLikes = zeros(parsed_args["nChroms"]+1)
				ref = nonref = 0
				nonref_allele ='N'

				ref = ((1:4)[mySite.reference .== alleles])[1]
				if ref == major
					nonref = minor
					nonref_allele = alleles[nonref]
				elseif ref == minor
					nonref = major
					nonref_allele = alleles[nonref]
				end

				# if nonref is defined and user wants to compute this quantities
				if nonref != 0 && parsed_args["nChroms"]>0

					for i = 0:parsed_args["nChroms"]
						safLikes[i+1] = calcFreqLogLike1_MajorMinor(myReads, ref, nonref, i/parsed_args["nChroms"], parsed_args["phredscale"])
					end
					safLikes = safLikes - maximum(safLikes)

					# calculate maximum (and expectation)
					freqMax = [((1:length(safLikes))[safLikes .== 0])]-1

					safLikes_norm = e.^safLikes
					safLikes_norm = safLikes_norm / sum(safLikes_norm)
					freqE = sum( (0:parsed_args["nChroms"]) .* safLikes_norm )
					safLikes_norm = 0

				end

				# write to files
				if parsed_args["fsaf"]!="/dev/null"
					write(f3, join(safLikes,"\t"), "\n") 
				end

				write(f2, join( (mySite.chrom, mySite.position, mySite.reference, nonref_allele, alleles[major], alleles[minor], lrtSNP, lrtBia, lrtTria, freqsMLE[2], freqMax/parsed_args["nChroms"], freqE/parsed_args["nChroms"]), "\t"), "\n")
				
				if parsed_args["verbose"]>0 && (allSites % parsed_args["printSites"])==0
					println("allSites:", allSites, "; ", "infSites:", infSites)
				end

			end
		end
	end

GZip.close(f2)

if parsed_args["fsaf"]!="/dev/null"
	GZip.close(f3)
end

if parsed_args["verbose"] > 0
	println("Done!\nallSites:", allSites, "; ", "infSites:", infSites)
end

end





