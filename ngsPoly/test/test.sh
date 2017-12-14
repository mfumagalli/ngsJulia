
POLY=~/Software/ngsPoly

# simulate 1 sample with ploidy 

NS=10000

Rscript $POLY/writePars.R -k 1 -n 10000 -p 1 > Data/test.pars

for MD in 2 5 10 20 50; # this is per haploid
do

	echo \# $MD

	> Data/test.$MD.mpileup
	offset=0
	for i in 1 2 3 4 5 5 4 3 2 1;
	do
		CD=$(($MD*$i))
		echo $CD
		Rscript $POLY/simulMpileup.R --copy ${i}x1 --sites $NS --depth $CD --qual 20 --ksfs 1 --ne 10000 --panc 1 --offset $offset >> Data/test.$MD.mpileup
		offset=$(($offset+$NS))
	done

	gzip Data/test.$MD.mpileup

	~/Software/julia/julia $POLY/ngsPoly.jl --fin Data/test.$MD.mpileup.gz --fpars Data/test.pars --fout Data/test.$MD.out.gz --nSamples 1 --ploidy 1-6 --debug 2 --verbose 0 --keepRef 1 | gzip > Data/test.$MD.genolikes.gz

done

for MD in 2 5 10 20 50;
do

	echo $MD

	Rscript calcPloidyWin.R Data/test.$MD.out.gz 1000 1000 0.01 0 NULL > Results/test.$MD.likes

	Rscript like2bay.R Results/test.$MD.likes > Results/test.$MD.priors

	Rscript calcPloidyWin.R Data/test.$MD.out.gz 1000 1000 0.01 0 Results/test.$MD.priors > Results/test.$MD.pprobs

done


