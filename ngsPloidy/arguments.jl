
# Specific functions for ngsPoly

import GZip
using ArgParse
using Combinatorics

# http://argparsejl.readthedocs.io/en/latest/argparse.html
# parsing for ngsPoly
function parse_commandline_poly()
	s = ArgParseSettings()

	@add_arg_table s begin

		"--fin"
			help = "input file gzipped mpileup"
			arg_type = AbstractString
			required = true

		"--fpars"
			help = "pars (pANC and geno priors) file"
			arg_type = AbstractString
			default = "NULL"

		"--fout"
			help = "output file gzipped text"
			arg_type = AbstractString
			default = "NULL"

		"--fglikes"
                        help = "genotype likelihoods file gzipped text"
                        arg_type = AbstractString
                        default = "NULL"

		"--nSamples"
			help = "number of samples"
			arg_type = Int64
			required = true

		"--ploidy"
			help = "ploidies to be tested, e.g. 2-5 or 1,3-5"
			arg_type = AbstractString
			required = false
			default = "1-8"

		"--keepRef"
			help = "keep reference as one possible allele"
			arg_type = Int64
			required = false
			default = 0

		"--callGeno"
			help = "call genotypes (highest likelihood/posterior probability)"
			arg_type = Int64
			required = false
			default = 0

		"--unif"
			help = "use a uniform prior for genotype probabilities"
			arg_type = Int64
			required = false
			default = 0

		"--thSnp"
			help = "chisquare for SNP calling"
			arg_type = Float64
			default = -Inf

		"--thBia"
			help = "chisquare for biallelic calling"
			arg_type = Float64
			default = -Inf

		"--thTria"
			help = "chisquare for triallelic (non) calling"
			arg_type = Float64
			default = Inf

		"--minMaf"
			help = "minimum allele frequency for SNP calling"
			arg_type = Float64
			default = -Inf

		"--minQ"
			help = "minimum base quality in phredscore"
			arg_type = Int64
			default = 13

		"--minNonMajorCount"
			help = "minimum non major base count"
			arg_type = Int64
			default = 0

                "--minNonMajorProportion"
			help = "minimum non major base proportion"
			arg_type = Float64
			default = 0.0

		"--minGlobalDepth"
			help = "minimum global depth"
			arg_type = Int64
			default = 1

		"--maxGlobalDepth"
			help = "maximum global depth"
			arg_type = Int64
			default = 10000

		"--minSampleDepth"
			help = "minimum sample depth"
			arg_type = Int64
			default = 0

		"--maxSampleDepth"
			help = "maximum sample depth"
			arg_type = Int64
			default = 10000

		"--minSamples"
			help = "minimum number of valid samples to retain site"
			arg_type = Int64
			default = 1

		"--nGrids"
			help = "grid density for grid-search estimation of allele frequencies"
			arg_type = Int64
			default = 0

		"--tol"
			help = "tolerance for GSS estimation of allele frequencies"
			arg_type = Float64
			default = 1e-5

		"--phredscale"
			help = "phredscale"
			arg_type = Int64
			default = 33

		"--verbose"
			help = "verbosity level"
			arg_type = Int64
			default = 0

    		"--debug"
			help = "debug for 1 sample"
			arg_type = Int64
			default = 0

		"--printSites"
			help =  "print on stdout every --printSites sites"
			arg_type = Int64
			default = 10000

	end

	return parse_args(s)

end


function parse_commandline_polyLite()

        s = ArgParseSettings()

        @add_arg_table s begin

                "--fin"
                        help = "input file gzipped mpileup"
                        arg_type = AbstractString
                        required = true

                "--fout"
                        help = "output file gzipped text"
                        arg_type = AbstractString
                        default = "NULL"

                "--fglikes"
                        help = "genotype likelihoods file gzipped text"
                        arg_type = AbstractString
                        default = "NULL"

                "--nSamples"
                        help = "number of samples"
                        arg_type = Int64
                        required = true

                "--ploidy"
                        help = "ploidies to be tested, e.g. 2-5 or 1,3-5"
                        arg_type = AbstractString
                        required = false
                        default = "1-8"

		"--keepRef"
                        help = "keep reference as one possible allele"
                        arg_type = Int64
                        required = false
                        default = 0

                "--thSnp"
                        help = "chisquare for SNP calling"
                        arg_type = Float64
                        default = -Inf

                "--thBia"
                        help = "chisquare for biallelic calling"
                        arg_type = Float64
                        default = -Inf

                "--thTria"
                        help = "chisquare for triallelic (non) calling"
                        arg_type = Float64
                        default = Inf

                "--minMaf"
                        help = "minimum allele frequency for SNP calling"
                        arg_type = Float64
                        default = 0.0

                "--minQ"
                        help = "minimum base quality in phredscore"
                        arg_type = Int64
                        default = 13

                "--minNonMajorCount"
                        help = "minimum non major base count"
                        arg_type = Int64
                        default = 0

                "--minNonMajorProportion"
                        help = "minimum non major base proportion"
                        arg_type = Float64
                        default = 0.0

                "--minGlobalDepth"
                        help = "minimum global depth"
                        arg_type = Int64
                        default = 1

                "--maxGlobalDepth"
                        help = "maximum global depth"
                        arg_type = Int64
                        default = 10000

                "--minSampleDepth"
                        help = "minimum sample depth"
                        arg_type = Int64
                        default = 0

                "--maxSampleDepth"
                        help = "maximum sample depth"
                        arg_type = Int64
                        default = 10000

                "--minSamples"
                        help = "minimum number of valid samples to retain site"
                        arg_type = Int64
                        default = 1

		"--nGrids"
                        help = "grid density for grid-search estimation of allele frequencies"
                        arg_type = Int64
                        default = 0

                "--tol"
                        help = "tolerance for GSS estimation of allele frequencies"
                        arg_type = Float64
                        default = 1e-3

                "--phredscale"
                        help = "phredscale"
                        arg_type = Int64
                        default = 33

                "--verbose"
                        help = "verbosity level"
                        arg_type = Int64
                        default = 0

        end

        return parse_args(s)

end




# parse cases (ploidies to be tested) from input
function parsePloidy(cases::AbstractString)

	a = split(cases, ",")
	b = []
	for i in 1:length(a)
		if length(a[i]) > 1
			b = [b; parse(Int, split(a[i], "-")[1]):parse(Int, split(a[i], "-")[2])]
		else
			b = [b; parse(Int,a[i])]
		end
	end

	# to go from Any to Int64
	c = zeros(Int64, length(b))
	for i in 1:length(c)
		c[i] = b[i]
	end


	return sort(c)

end

# read anc/der and genotype probabilities in input
function readPriors(fname::AbstractString, ploidy::Array{Int64,1})

	# genotype priors
	pGenos = zeros(Float64, length(ploidy), length(ploidy)+1)
	# probs that major is ancestral or derived
	pMajorAnc = zeros(Float64, 2)

	open(fname) do file
		iter = -1
		for line in eachline(file)
			iter += 1
			if (iter <= length(ploidy))
				l = split(chomp(line), "\t")
				if iter==0
					# first line is probs that major is ancestral or derived
					for i in 1:2
						pMajorAnc[i] = parse(Float64, l[i])
					end
					pMajorAnc[1:2] = pMajorAnc[1:2] / sum(pMajorAnc[1:2])
				else
					for i in 1:(iter+1)
						pGenos[iter,i] = parse(Float64, l[i])
					end
					pGenos[iter,1:(length(ploidy)+1)] = pGenos[iter,1:(length(ploidy)+1)] / sum(pGenos[iter,1:(length(ploidy)+1)])
				end
			end
		end
	end

	return (pMajorAnc, pGenos)
end

# calculate AIC and BIC from model likelihoods # deprecated, not working
function calcModelStat(likes::Array{Float64,2}, ploidy::Array{Int64, 1}, sites::Array{Float64,1})

	if length(ploidy)!=size(likes,2)
		println("Error!: ploidy and likelihoods size do not match in AIC/BIC function.")
	end

	nsams = size(likes, 1)

	# nr of parameters, set of genotypes
	k = ploidy +1
	# fix this in case SNP calling is done (k is less) or it is unfolded (k = ploidy) or ancestral uncertain (k=...)

	bic = zeros(Float64, nsams, length(ploidy))
	aic = zeros(Float64, nsams, length(ploidy))

	for n = 1:nsams

		for p = 1:length(ploidy)

			# this gives error, no method matching +(::Vector{Int64}, ::Int64)
			bic[n,p] = -2*likes[n,p]+(k[p]*log(sites[n]))
			aic[n,p] = (2*k[p])-(2*likes[n,p]) + ((2*k[p]*(k[p]+1))/(sites[n]-k[p]-1))

		end
	end

	return (aic, bic)

end

# ---------------------------------------------------
