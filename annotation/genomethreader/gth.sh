#!/bin/bash
#SBATCH --job-name=gth
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=basic
#SBATCH --output=slurm.%x.out
#SBATCH --nice=5000
#SBATCH --time=1-00:00:00

module load conda
conda activate annotation

GENOME=${1:-Cmultiflora.genome.fna.gz}
PROT_PREFIX=${2:-proteins}
OUT_PREFIX=${3:-Cmultiflora.gth}
PART=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})

PROTEINS=$(ls "$PROT_PREFIX".part_"$PART".faa*) #proteins.part_001.faa.gz
OUT="$OUT_PREFIX".part_"$PART".gff3 #Cmultiflora.gth.part_001.out

rsync --copy-links --progress "$GENOME" "$PROTEINS" "$TMPDIR"
pushd "$TMPDIR"

stdbuf -i0 -o0 -e0 gth -genomic "$GENOME" -protein "$PROTEINS" -gff3out -intermediate -v -o "$OUT"

popd
rsync --progress "$TMPDIR"/"$OUT" .
