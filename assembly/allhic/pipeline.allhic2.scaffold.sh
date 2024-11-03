#!/bin/bash
#SBATCH --job-name=allhic-optimize
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --partition=basic
#SBATCH --constraint=array-1core|array-4core
#SBATCH --nice=500
#SBATCH --output=slurm.%x.out
#SBATCH --time=2-00:00:00

export PATH=/home/user/kramml/git/ALLHiC/scripts/:/home/user/kramml/git/ALLHiC/bin/:$PATH
module load samtools

PREFIX=$1
K=${2:=30}

allhic optimize $PREFIX.counts_GATC."$K"g${SLURM_ARRAY_TASK_ID}.txt $PREFIX.clm
