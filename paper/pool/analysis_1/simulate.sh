

#SIM=/home/mfumagal/Software/HMMploidy

SIM=.

NSITES=100000

for S in 10 20 50 100;
do
	for D in 0.1 0.5 1 2 5 10;
	do

		echo $S $D

		Rscript $SIM/simulMpileup.R --out Data/$S.$D.txt --copy 2x$S --sites $NSITES --depth $D --errdepth 0 --qual 20 --ksfs 1 --ne 10000 --pool | gzip > Data/$S.$D.mpileup.gz

	done
done







