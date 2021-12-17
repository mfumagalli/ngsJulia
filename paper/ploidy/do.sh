
JULIA=~/Software/julia-1.6.1/bin/julia
NGSJULIA=~/Software/ngsJulia

Rscript $NGSJULIA/simulMpileup.R --copy 2x1,3x8,4x1 --sites 1000 --depth 10 --ne 10000 --k 1 | gzip > test.mpileup.gz

Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 1 > test.pars
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p 0.5 > test.fold.pars
Rscript $NGSJULIA/ngsPloidy/writePars.R -k 1 -n 10000 -p -1 > test.auto.pars

echo 1
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.mpileup.gz --fpars test.pars --nSamples 10 --keepRef 1 > unfold.txt

echo 2
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.mpileup.gz --fpars test.auto.pars --nSamples 10 --keepRef 1 > auto.txt

echo 2
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.mpileup.gz --fpars test.fold.pars --nSamples 10 --keepRef 1 > fold.txt

echo 3
$JULIA $NGSJULIA/ngsPloidy/ngsPloidy.jl --fin test.mpileup.gz --nSamples 10 > sample.txt


