#!/bin/bash
#SBATCH --job-name=pasa
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=slurm.%x.out
#SBATCH --time=1-00:00:00

module load perl
export PATH=/home/user/kramml/git/PASApipeline-v2.5.2/:/home/user/kramml/git/PASApipeline-v2.5.2/bin:$PATH

GENOME=${1:-Cmultiflora.genome.diploid.fna}
TRANSCRIPTS=${2:-Cmultiflora.transcripts.fna.clean}
CONFIG=${3:-alignAssembly.config}

Launch_PASA_pipeline.pl -c "$CONFIG" --ALT_SPLICE \
                        -g "$GENOME" \
                        -t "$TRANSCRIPTS" \
                        --CPU 4
