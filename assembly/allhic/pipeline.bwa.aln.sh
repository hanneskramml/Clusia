#!/bin/bash
#
#SBATCH --job-name=bwa-aln
#SBATCH --cpus-per-task=24
#SBATCH --mem=30G
#SBATCH --mail-type=end
#SBATCH --output=slurm.%x.%j.out
#SBATCH --time=2-00:00:00

module load bwa samtools

GENOME=$1
READS=$2
OUT=$3 # e.g. *.alignment.R1.sai

# bwa index -a bwtsw $GENOME
# samtools faidx $GENOME

rsync --copy-links --progress $GENOME $GENOME.* $READS $TMPDIR/
bwa aln -t 24 $TMPDIR/$GENOME $TMPDIR/$READS  > $OUT
