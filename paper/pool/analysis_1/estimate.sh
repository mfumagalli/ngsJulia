
JULIA=/home/mfumagal/Software/ngsJulia

POOL=$JULIA/ngsPool

J=~/Software/julia-1.6.0/bin/julia

for S in 10 20 50 100;
do
	for D in 0.1 0.5 1 2 5 10;
	do

		echo $S $D

		$J $POOL/ngsPool.jl --fin Data/$S.$D.mpileup.gz --fout Results/$S.$D.out.gz --nSamp `expr $S + $S` --fsaf Results/$S.$D.saf.gz --lrtSnp -Inf >> /dev/null

	done
done

	
# then pick one case to do QUAL=15
# then pick one case to increase Ne or play with ksfs






