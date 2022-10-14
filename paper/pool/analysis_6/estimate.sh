
JULIA=/home/mfumagal/Software/ngsJulia

POOL=$JULIA/ngsPool

J=~/Software/julia-1.6.0/bin/julia

NSITES=100

SS=150

for S in cases controls;
do
	for F in null causal;
	do
		for D in 0.5 1 2 5;
		do

			echo $S $F $D

			$J $POOL/ngsPool.jl --fin Data/$S.$F.$D.mpileup.gz --fout Results/$S.$F.$D.out.gz --nSamp `expr $SS + $SS` --fsaf Results/$S.$F.$D.saf.gz --lrtSnp -Inf >> /dev/null

		done
	done
done







