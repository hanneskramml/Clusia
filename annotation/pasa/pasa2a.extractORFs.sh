#!/bin/bash
#SBATCH --job-name=pasa
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=slurm.%x.out
#SBATCH --time=12:00:00

module load conda perl samtools
export PATH=/home/user/kramml/git/PASApipeline-v2.5.2/:/home/user/kramml/git/PASApipeline-v2.5.2/bin:$PATH

PASA_FASTA=${1:-clusia.assemblies.fasta}
PASA_GFF=${2:-clusia.pasa_assemblies.gff3}

# Extraction of ORFs from PASA assemblies
./pasa/scripts/pasa_asmbls_to_training_set.dbi \
       --pasa_transcripts_fasta "$PASA_FASTA" \
       --pasa_transcripts_gff3 "$PASA_GFF" \
       -S
