#!/bin/bash
#
#SBATCH --job-name=minimap2
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
#SBATCH --mail-type=end
#SBATCH --output=slurm.%x.out
#SBATCH --time=0-06:00:00

module load conda samtools
conda activate hic

GENOME=$1
READS=$2
MODE=${3:-map-pb} # Map long noisy genomic reads: map-pb | Full genome/assembly alignment: asm5/asm10/asm20 - for ~0.1/1/5% sequence divergence. "asm5" for sequence divergence below 1%, "asm10" for divergence around a couple of percent and "asm20" for divergence not more than 10%
OUT=${4:-${GENOME}.alignment.bam}

rsync --copy-links --progress "$GENOME" "$READS" "$TMPDIR"/
#minimap2 -ax map-pb -t 16 $TMPDIR/$GENOME $TMPDIR/$READS --secondary=no | samtools sort -@ 8 -m 1G -o $GENOME.alignment.bam -T $TMPDIR/tmp.ali   # A secondary alignment occurs when a given read could align reasonably well to more than one place. purge_haplotigs
minimap2 -ax "$MODE" -t 16 "$TMPDIR"/"$GENOME" "$TMPDIR"/"$READS" | samtools sort -@ 8 -m 1G -T "$TMPDIR"/tmp.ali -o "$OUT"
