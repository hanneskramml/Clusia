#!/bin/bash
#SBATCH --job-name=pasa
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=slurm.%x.transdb.out
#SBATCH --time=10:00:00

module load conda perl samtools
export PATH=/home/user/kramml/git/PASApipeline-v2.5.2/:/home/user/kramml/git/PASApipeline-v2.5.2/bin:$PATH

TRANSCRIPTS=${1:-Cmultiflora.transcripts.fna.clean}
CONFIG=${2:-alignAssembly.config}

# generate the comprehensive transcriptome database
./pasa/scripts/build_comprehensive_transcriptome.dbi \
           -c "$CONFIG" \
           -t "$TRANSCRIPTS" \
           --min_per_ID 95 \
           --min_per_aligned 30
