 #define the type and hierarchy

# struct Read # commented out as not used in generic.jl (2020/5/16)
# 	base::Char
# 	baseQuality::Char
# end
# #Read =
# Read('N', '!')

mutable struct Reads #mutable added
        base::AbstractString
        baseQuality::AbstractString
end
# Reads =
Reads("N", "!")

mutable struct Site #mutable added(2020/4/17) for convertSyms(myReads, mySite) line36
	# chrom::AbstractString
	chrom::SubString #as how to convert SubString{String} to AbstractString
	position::Int64
	reference::Char
end
# Site =
Site("chrom", 0, 'N')

# --------------------

abstract type Sequence end

struct Genotype <: Sequence
	ploidy::Int64
	allele::AbstractString
	logLike::Float64
end
#Genotype =
Genotype(0, "N", -1.0)

struct Frequency <: Sequence # Frequency is a subtype of Sequence
	frequency::Float64
	probability::Float64
end
# Frequency =
Frequency(0,0)

# --------------------------
