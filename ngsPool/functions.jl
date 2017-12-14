
# Specific Functions for ngsPool

import GZip
using ArgParse

# http://argparsejl.readthedocs.io/en/latest/argparse.html
function parse_commandline_pool()
	s = ArgParseSettings()

	@add_arg_table s begin

		"--fin"
			help = "input file gzipped mpileup"
			arg_type = AbstractString
			required = true

		"--fout"
			help = "output file gzipped text"
			arg_type = AbstractString
			required = true

		"--fsaf"
			help = "output gzipped saf file"
			arg_type = AbstractString
			default = "/dev/null"

		"--nChroms"
			help = "number of samples [>0 ensables saf likelihoods]"
			arg_type = Int64
			default = 0

		"--thSnp"
			help = "chisquare for SNP calling"
			arg_type = Float64
			default = 7.82

		"--thBia"       
			help = "chisquare for biallelic calling"
			arg_type = Float64
			default = -Inf
		
		"--thTria"
			help = "chisquare for triallelic (non) calling"
			arg_type = Float64
			default = Inf

		"--minQ"
			help = "minimum base quality in phredscore"
			arg_type = Int64
			default = 5

		"--minDepth"
			help = "minimum global depth"
			arg_type = Int64
			default = 1

		"--maxDepth"
			help = "maximum global depth"
			arg_type = Int64
			default = 100000

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
			default = 1

		"--printSites"
			help =  "print on stdout every --printSites sites"
			arg_type = Int64
			default = 10000

	end

	return parse_args(s)

end

# -----------------------------------------

