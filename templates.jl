
global ALLELES = ['A', 'C', 'G', 'T'];

# define the type and hierarchy

mutable struct Reads #mutable added
        base::AbstractString
        baseQuality::AbstractString
end

Reads("N", "!");

mutable struct Site
	chrom::SubString #as how to convert SubString{String} to AbstractString
	position::Int64
	reference::Char
end

Site("chrom", 0, 'N');

# --------------------

abstract type Sequence end

struct Genotype <: Sequence
	ploidy::Int64
	allele::AbstractString
	logLike::Float64
end

Genotype(0, "N", -1.0);

struct Frequency <: Sequence # Frequency is a subtype of Sequence
	frequency::Float64
	probability::Float64
end

Frequency(0,0);

# --------------------------



