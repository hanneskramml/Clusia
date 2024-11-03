#!/bin/bash
#
#SBATCH --job-name=repmask2
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --output=slurm.%x.%j.out
#SBATCH --time=8-00:00:00

# barely any plant TEs covered in DFAM DB => use ensembl DB instead (nrTEplants)
# nrTEplants DB takes about a week to complete for C. multiflora, possibly submit every chromosome separately or bad config:
# -pa(rallel) number of instances. rmblast uses 4 threads per instance.

module load repeatmasker	#/4.1.3-p1
conda activate repeatmasker

GENOME=$1
LIB=$2
OUT_PREFIX=$3

rsync --copy-links --progress "$GENOME" "$LIB" "$TMPDIR"
pushd "$TMPDIR"

#RepeatMasker -e rmblast -pa 20 -xsmall -a -nolow -species "$DB" data/"$GENOME"
RepeatMasker -e rmblast -pa 7 -xsmall -a -nolow -lib "$LIB" "$GENOME"
rename -v "$GENOME" "$OUT_PREFIX" "$GENOME".*
rename -v .masked .fna "$OUT_PREFIX"*

popd
rsync --progress "$TMPDIR"/"$OUT_PREFIX"* .
