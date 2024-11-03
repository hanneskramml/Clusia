#!/bin/bash
#SBATCH --job-name=stringtie
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=basic
#SBATCH --nice=0
#SBATCH --output=slurm.%x.out
#SBATCH --time=1-00:00:00

module load stringtie

ALN=${1:-Cmultiflora.alignment.bam}
OUT=${2:-Cmultiflora.transcripts.gtf}
LABEL=${3:-Cmu.STRG}

stringtie --rf -p 2 -l "$LABEL" -o "$OUT" "$ALN"
