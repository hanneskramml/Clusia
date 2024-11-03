#!/bin/bash
#SBATCH --job-name=liftoff
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=basic
#SBATCH --output=slurm.%x-%j.out
#SBATCH --time=04:00:00

# AllelicGenes: 1 CPU, 10G mem, 1h
# CmiCro2Cmu:   1 CPU, 15G mem, 4h

module load conda
conda activate annotation

GFF=${1:-Cmultiflora_v2.2.annotation.allelicGenes.gff3.gz}
REF=${2:-Cmultiflora_v2.haplotigs.softmasked.fna.gz}
TARGET=${3:-Cmultiflora_v2.scaffolds.softmasked.fna.gz}
OUT_PREFIX=${4:-Cmultiflora_v2.2.annotation.allelicGenes.chrom_choords}

if [ "${REF##*.}" == "gz" ]; then
   if [ ! -f "${REF%.gz}" ]; then pigz -vdfkp1 "$REF"; fi
   REF="${REF%.gz}"
fi

if [ "${TARGET##*.}" == "gz" ]; then
   if [ ! -f "${TARGET%.gz}" ]; then pigz -vdfkp1 "$TARGET"; fi
   TARGET="${TARGET%.gz}"
fi

liftoff -g "$GFF" -o "$OUT_PREFIX".gff3 -u "$OUT_PREFIX".unmapped_features.txt -dir "$OUT_PREFIX".intermediate_files -copies -sc 0.98 -flank 0.1 -polish "$TARGET" "$REF"
