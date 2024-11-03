#!/bin/bash
#
#SBATCH --job-name=orthofinder
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm.%x-%j.out

# eudicots: 16 CPUs, 64G, 2D
# malpighiales: 8 CPUs, 40G, 1D
# clusia: 8 CPUs, 20G, 10H

module load orthofinder mafft conda
conda activate fasttree-2.1.11

DATASET=${1:-eudicots}
INPUT="$DATASET"/dataset
OUTPUT="$TMPDIR"/output

echo "Running OrthoFinder... (dataset: "$DATASET")"
orthofinder.py -S diamond -M msa -A mafft -T fasttree -X -t 16 -p "$TMPDIR" -f "$INPUT" -o "$OUTPUT"

echo "Archiving/cleaning files..."
./cleanup.sh "$OUTPUT"/Results_*

echo "Moving final output..."
rsync -rt --info=progress2 "$OUTPUT"/Results_*/* "$DATASET"
