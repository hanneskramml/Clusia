#!/bin/bash
#
#SBATCH --job-name=kmc
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=end
#SBATCH --output=slurm.%x.out

echo "$TMPDIR/Crosea.poly.trimmedreads.fna.gz" > $TMPDIR/reads.fofn
echo "$TMPDIR/Crosea.preads.fna.gz" >> $TMPDIR/reads.fofn
cat $TMPDIR/reads.fofn

rsync --copy-links --progress Crosea.preads.fna.gz Crosea.poly.trimmedreads.fna.gz $TMPDIR/

# kmer 21, 16 threads (only utilizes one CPU?, no!), 64G of memory, counting kmer coverages between 1 and 1e6x (to cover repetitive regions)
kmc/kmc -k21 -t16 -m64 -ci1 -cs1000000 -fm @$TMPDIR/reads.fofn Crosea.mergedreads.kmcdb $TMPDIR
