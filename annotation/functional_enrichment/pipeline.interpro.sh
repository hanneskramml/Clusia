#!/bin/bash
#
#SBATCH --job-name=interpro
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
#SBATCH --output=slurm.%x.out
#SBATCH --time=1-00:00:00
#SBATCH --license=interpro

QUERY=${1:-Cmultiflora.cds.faa}
OUT=${2:-Cmultiflora.interpro.tsv}

/scratch/mirror/interpro/interproscan-5.62-94.0/interproscan.sh -i "$QUERY" -T "$TMPDIR" -t p --goterms --iprlookup --cpu 16 --formats TSV --outfile "$OUT"
