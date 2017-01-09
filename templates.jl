
# Templates

type Read
	base::Char
	baseQuality::Char
end
Read() = Read('N', '!');

type Reads
        base::AbstractString
        baseQuality::AbstractString
end
Reads() = Reads("N", "!");

type Site
	chrom::AbstractString
	position::Int64
	reference::Char
end
Site() = Site("chrom", 0, 'N');

# --------------------

abstract Sequence

type Genotype <: Sequence
	ploidy::Int64
	allele::AbstractString
	logLike::Float64
end
Genotype() = Genotype(0, "N", -1.0);

type Frequency <: Sequence # Frequency is a subtype of Sequence
	frequency::Float64
	probability::Float64
end
Frequency() = Frequency(0,0)

# --------------------------


