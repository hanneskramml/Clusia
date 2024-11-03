#!/bin/bash
#
#SBATCH --job-name=eggnog
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=slurm.%x.out
#SBATCH --time=02:00:00

module load conda
conda activate eggnog-mapper-2.1.10

export EGGNOG_DATA_DIR=./db     # database not globally available, download with: download_eggnog_data.py --data_dir db

PROT=${1:-Cmultiflora.cds.faa}
OUT_PREFIX=${2:-Cmultiflora}

# 2-step mode (search + annotation). Search freates $OUT_PREFIX.emapper.seed_orthologs
emapper.py -i "$PROT" -m diamond --dmnd_iterate yes --sensmode ultra-sensitive --no_annot --cpu 8 --temp_dir "$TMPDIR" --output "$OUT_PREFIX"
