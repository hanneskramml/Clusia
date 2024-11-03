#!/bin/bash
#
#SBATCH --job-name=repmask3
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --output=slurm.%x.%j.out
#SBATCH --time=0-12:00:00

module load repeatmasker	#/4.1.3-p1
conda activate repeatmasker

GENOME=$1
LIB=$2
OUT_PREFIX=$3

rsync --copy-links --progress "$GENOME" "$LIB" "$TMPDIR"
pushd "$TMPDIR"

RepeatMasker -e rmblast -pa 10 -xsmall -a -nolow -lib "$LIB" "$GENOME"
rename -v "$GENOME" "$OUT_PREFIX" "$GENOME".*
rename -v .masked .fna "$OUT_PREFIX"*

popd
rsync --progress "$TMPDIR"/"$OUT_PREFIX"* .
