#!/bin/bash
#
#SBATCH --job-name=busco
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --output=slurm.%x.%j.out
#SBATCH --time=0-12:00:00

module load conda
conda activate busco

IN=$1
OUT=${2:-${1%.fna*}.busco}

busco -i "$IN" -m genome -l eudicots_odb10 --download_path "$HOME/.busco" -c4 -o "$OUT"
mv "$OUT"/short_summary.*.txt "$OUT".txt
