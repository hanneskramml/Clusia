#!/bin/bash
#SBATCH --job-name=pasa
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --output=slurm.%x.out
#SBATCH --time=2-00:00:00

module load perl
export PATH=/home/user/kramml/git/PASApipeline-v2.5.2/:/home/user/kramml/git/PASApipeline-v2.5.2/bin:/home/user/kramml/git/fasta-36.3.8i/bin:$PATH

GFF=${1:-Cmultiflora.annot0.gff3}
GENOME=${2:-Cmultiflora.genome.diploid.fna}
TRANSCRIPTS=${3:-Cmultiflora.transcripts.fna.clean}
CONFIG=${4:-annotCompare.config}

# validate gene models
./pasa/misc_utilities/pasa_gff3_validator.pl "$GFF"

# Load annotations
# ./pasa/scripts/Load_Current_Gene_Annotations.dbi \
#      -c "$CONFIG" -g "$GENOME" -P "$GFF"

# Annotation comparision - run at least two cycles
Launch_PASA_pipeline.pl \
        -c "$CONFIG" -A \
        -g "$GENOME" -t "$TRANSCRIPTS" \
        -L --annots "$GFF" \
        --CPU 8
