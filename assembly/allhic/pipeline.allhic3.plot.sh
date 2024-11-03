#!/bin/bash
#SBATCH --job-name=allhic-plot
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=basic
#SBATCH --constraint=array-1core
#SBATCH --output=slurm.%x.out
#SBATCH --time=1-00:00:00

export PATH=/home/user/kramml/git/ALLHiC/scripts/:/home/user/kramml/git/ALLHiC/bin/:$PATH

PREFIX=$1
K=${2:-30}

mkdir "$PREFIX".scaffolds.k"$K".plots
cd ./"$PREFIX".scaffolds.k"$K".plots
ALLHiC_plot ../"$PREFIX".bam ../"$PREFIX".scaffolds.k"$K".agp ../"$PREFIX".scaffolds.k"$K".chrom.sizes 500k pdf
