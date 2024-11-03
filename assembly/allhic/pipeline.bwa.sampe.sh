#!/bin/bash
#
#SBATCH --job-name=bwa-sampe
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --mail-type=end
#SBATCH --output=slurm.%x.%j.out
#SBATCH --time=2-00:00:00

module load bwa
GENOME=$1
ALN1=$2
ALN2=$3
READS1=$4
READS2=$5
OUT=$6

rsync --copy-links --progress $GENOME $GENOME.* $ALN1 $ALN2 $READS1 $READS2 $TMPDIR/
bwa sampe $TMPDIR/$GENOME $TMPDIR/$ALN1 $TMPDIR/$ALN2 $TMPDIR/$READS1 $TMPDIR/$READS2 > $OUT
