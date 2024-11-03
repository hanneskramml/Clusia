#!/bin/bash
#
#SBATCH --job-name=genomescope
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --mail-type=end
#SBATCH --output=slurm.%x.out

PREFIX=$1 #Cmultiflora.hicreads.kmcdb
echo $PREFIX

# sbatch slurm.kmc.sh $READS
kmc/kmc_tools transform $PREFIX histogram $PREFIX.hist_k21_cx1e6 -cx1000000

# conda activate hic
#L=$(smudgeplot.py cutoff $PREFIX.hist_k21_cx1e6 L) #  20 -  200
L=20
U=$(smudgeplot.py cutoff $PREFIX.hist_k21_cx1e6 U) # 500 - 3000
echo $L #Cro: 10
echo $U #Cro: 1300

kmc/kmc_tools transform $PREFIX -ci"$L" -cx"$U" reduce $PREFIX.L"$L"_U"$U"
kmc/smudge_pairs $PREFIX.L"$L"_U"$U" $PREFIX.L"$L"_U"$U".coverages.tsv $PREFIX.L"$L"_U"$U".pairs.tsv > $PREFIX.L"$L"_U"$U".familysizes.tsv

# without smudge_pairs
#kmc/kmc_tools transform Cmultiflora.preads.fna.gz.kmcdb -ci10 -cx1300 dump -s Cmultiflora.kmcdb_L10_U1300.dump
#srun -c1 --mem=200G smudgeplot.py hetkmers -o Cmultiflora.kmcdb_L10_U1300 < Cmultiflora.kmcdb_L10_U1300.dump

# smudgeplot
smudgeplot.py plot -t $PREFIX -k21 -o $PREFIX.L"$L"_U"$U". $PREFIX.L"$L"_U"$U".coverages.tsv

# genomescope, adjust ploidy lever (-p)
Rscript genomescope.R -i $PREFIX.hist_k21_cx1e6 -n $PREFIX.genscope. -k21 -p4 -o .

#Rscript genomescope.R -i Cmultiflora.poly.trimmedreads.kmcdb.hist_k21_cx1e6 -n Cmultiflora.trimmedreads.genscope. -k21 -p2 -o .
#Rscript genomescope.R -i Cmultiflora.poly.trimmedreads.kmcdb.hist_k21_cx1e6 -n Cmultiflora.trimmedreads.genscope.p4 -k21 -p4 -l47 -o .
