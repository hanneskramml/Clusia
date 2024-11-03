#!/bin/bash
#
#SBATCH --job-name=purgehaplotigs3
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --mail-type=end
#SBATCH --output=slurm.%x-%j.out

# conda activate hic
SPECIES=$1

#rsync --copy-links --progress $SPECIES.assembly.fna $SPECIES.assembly.fna.fai $SPECIES.alignment.bam $SPECIES.alignment.bam.bai $TMPDIR/
purge_haplotigs purge -g $SPECIES.assembly.fna -c $SPECIES.alignment.bam.coverage_stats.csv -r $SPECIES.repeats.bed -d -b $SPECIES.alignment.bam -t 8 -f -o $SPECIES.curated
