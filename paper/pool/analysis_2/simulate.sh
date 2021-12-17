

NSITES=10000

SIM=.

for S in 20 50;
do
	for F in 0.1 0.2 0.3 0.4 0.5;
	do
		for D in 0.5 1 2 5;
		do

			echo $S $F $D

			Rscript $SIM/simulMpileup_qq.R --out Data/$S.$F.$D.txt --copy 2x$S --sites $NSITES --depth $D --errdepth 0 --qual 20 --ksfs 1 --ne 10000 --qq $F --pool | gzip > Data/$S.$F.$D.mpileup.gz

		done
	done
done







