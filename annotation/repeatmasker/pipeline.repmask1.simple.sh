#!/bin/bash
#
#SBATCH --job-name=repmask1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --output=slurm.%x.%j.out
#SBATCH --time=0-03:00:00

module load repeatmasker	#/4.1.3-p1 Cmu.scaffolds
conda activate repeatmasker

GENOME=$1
OUT_PREFIX=$2

rsync --copy-links --progress "$GENOME" "$LIB" "$TMPDIR"
pushd "$TMPDIR"

RepeatMasker -e rmblast -pa 8 -xsmall -a -noint "$GENOME"
rename -v "$GENOME" "$OUT_PREFIX" "$GENOME".*
rename -v .masked .fna "$OUT_PREFIX"*

popd
rsync --progress "$TMPDIR"/"$OUT_PREFIX"* .
