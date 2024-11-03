#!/bin/bash
#
#SBATCH --job-name=purgehaplotigs1
#SBATCH --cpus-per-task=8
#SBATCH --mem=1G
#SBATCH --mail-type=end
#SBATCH --output=slurm.%x.%j.out
#SBATCH --time=1-00:00:00

# module unload perl
# conda activate hic
# samtools index alignment.bam

GENOME=$1
ALIGNMENT=$2
CUTOFF="${3:-400}" # cutoff for read depth

purge_haplotigs hist -b $ALIGNMENT -g $GENOME -d $CUTOFF -t 8
