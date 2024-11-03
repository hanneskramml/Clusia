#!/bin/bash
#
#SBATCH --job-name=orthofinder
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm.%x-%j.out

# eudicots: 8 CPUs, 20G, 3D
# clusia.vitis: 4 CPUs, 20G, 6H
# clusia.subgenomes: 4 CPUs, 20G, 4h

module load orthofinder

DATASET=${1:-eudicots}
INPUT="$DATASET"/dataset
OUTPUT="$TMPDIR"/output

echo "Running OrthoFinder... (dataset: "$DATASET")"
orthofinder.py -S diamond -X -t 8 -p "$TMPDIR" -f "$INPUT" -o "$OUTPUT"

echo "Archiving/cleaning files..."
./cleanup.sh "$OUTPUT"/Results_*

echo "Moving final output..."
rsync -rt --info=progress2 "$OUTPUT"/Results_*/* "$DATASET"
#find  . -name "*SpeciesIDs*" -exec cp {} $HOMEFOLDER/results \;
