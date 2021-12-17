
# Whensimulating causal SNPs in association studies, MAFs forcases and controls were assigned using a multiplicativedisease model. For this model, the prevalence of the dis-ease was fixed at 10%. We examined two sets of MAFsand relative risks. First, the combined MAF in cases andcontrols was 1% and the relative risk was 2. Second, thecombined MAF was 5% and the relative risk was 1.5. Asan example, with a combined MAF of 1% and a relativerisk of 2.0, the obtained MAFs for cases and controls are1.98% and 0.89%, respectively.

#http://zzz.bwh.harvard.edu/cgi-bin/cc2k.cgi

#Case-control parameters
#Number of cases 	150
#Number of controls 	150
#High risk allele frequency (A) 	0.1
#Prevalence 	0.2
#Genotypic relative risk Aa 	2
#Genotypic relative risk AA 	4
#Genotypic risk for aa (baseline) 	0.1653

#Marker locus B
#High risk allele frequency (B) 	0.05
#Penetrance at marker genotype bb 	0.1831
#Penetrance at marker genotype Bb 	0.348
#Penetrance at marker genotype BB 	0.6612
#Genotypic odds ratio Bb 	2.38
#Genotypic odds ratio BB 	8.703

#Expected allele frequencies
#	Case	Control
#B	0.09091	0.03977

#Case-control statistics: allelic 1 df test (B versus b)
#N cases for 80% power
#Alpha 	Power 	N cases for 80% power
#0.1 	0.8131 	144
#0.05 	0.7171	183 


SIM=.
NSITES=100

# cases
for S in 150;
do
	for F in 0.1;
	do
		for D in 0.5 1 2 5;
		do
			echo $S $F $D
			Rscript $SIM/simulMpileup_qq.R --out Data/cases.null.$D.txt --copy 2x$S --sites $NSITES --depth $D --errdepth 0 --qual 20 --ksfs 1 --ne 10000 --qq $F --pool | gzip > Data/cases.null.$D.mpileup.gz
		done
	done
done
# cases causual
for S in 150;
do
        for F in 0.09;
        do
                for D in 0.5 1 2 5;
                do
                        echo $S $F $D
                        Rscript $SIM/simulMpileup_qq.R --out Data/cases.causal.$D.txt --copy 2x$S --sites $NSITES --depth $D --errdepth 0 --qual 20 --ksfs 1 --ne 10000 --qq $F --pool | gzip > Data/cases.causal.$D.mpileup.gz
                done
        done
done

# controls
for S in 150;
do
        for F in 0.1;
        do
                for D in 0.5 1 2 5;
                do
                        echo $S $F $D
                        Rscript $SIM/simulMpileup_qq.R --out Data/controls.null.$D.txt --copy 2x$S --sites $NSITES --depth $D --errdepth 0 --qual 20 --ksfs 1 --ne 10000 --qq $F --pool | gzip > Data/controls.null.$D.mpileup.gz
                done
        done
done
# controls causual
for S in 150;
do
        for F in 0.04;
        do
                for D in 0.5 1 2 5;
                do
                        echo $S $F $D
                        Rscript $SIM/simulMpileup_qq.R --out Data/controls.causal.$D.txt --copy 2x$S --sites $NSITES --depth $D --errdepth 0 --qual 20 --ksfs 1 --ne 10000 --qq $F --pool | gzip > Data/controls.causal.$D.mpileup.gz
                done
        done
done



