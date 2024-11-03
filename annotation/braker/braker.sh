#!/bin/bash
#SBATCH --job-name=braker
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --mail-type=end
#SBATCH --output=slurm.%x.out
#SBATCH --time=10-00:00:00

# setup environment
module load genemark-es/4.65
conda activate braker2-2.1.6
export AUGUSTUS_CONFIG_PATH=/scratch/ecogenomics/clusia/braker/augustus
export PROTHINT_PATH=/apps/genemark-es/4.65/ProtHint/bin

# input parameters
GENOME=${1:-Cmultiflora.genome.fna}
SPECIES=${2:-Clusia_multiflora}
PROT=${3:-proteins.faa}
DIR=${4:-GENOME%.fna*}

# run braker pipeline C
braker.pl --species="$SPECIES" --genome="$GENOME" --prot_seq="$PROT" --cores=10 --workingdir="$DIR" --gff3 --softmasking
