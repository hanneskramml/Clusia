#!/bin/bash
#
#SBATCH --job-name=kmc
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=end
#SBATCH --output=slurm.%x.out

READS=$1
rsync --copy-links --progress $READS $TMPDIR/

# kmer 21, 16 threads, 64G of memory, counting kmer coverages between 1 and 1e6x (to cover repetitive regions)
kmc/kmc -k21 -t16 -m64 -ci1 -cs1000000 -fm $TMPDIR/$READS $READS.kmcdb $TMPDIR
