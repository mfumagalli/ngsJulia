
echo simulate
bash simulate.sh

echo estimate
bash estimate.sh

echo calculate
Rscript calculate.R > results.1.txt

#Rscript plot.1.R


