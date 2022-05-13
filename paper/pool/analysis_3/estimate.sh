
JULIA=/home/mfumagal/Software/ngsJulia

POOL=$JULIA/ngsPool

J=~/Software/julia-1.6.0/bin/julia

NSITES=10000

for S in 20 50;
do
	for F in 0 0.01 0.02 0.025 0.05;
	do
		for D in 0.5 1 2 5;
		do

			echo $S $F $D

			$J $POOL/ngsPool.jl --fin Data/$S.$F.$D.mpileup.gz --fout Results/$S.$F.$D.out.gz --nSamp `expr $S + $S` --fsaf Results/$S.$F.$D.saf.gz --lrtSnp -Inf >> /dev/null

		done
	done
done






