#!/bin/bash

#SBATCH --job-name=wgd.syn
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=0-00:10:00
#SBATCH --output=slurm.%x-%j.out

# input parameter
FAMILIES=${1:-wgd.dmd/Clusia_multiflora_H1.paranome.tsv}
KS=${2:-wgd.ksd/Clusia_multiflora_H1.paranome.tsv.ks.tsv}
GFF=${3:-Clusia_multiflora_H1.gff3.gz}
OUTDIR=${4:-wgd.syn.H1}

# load/set envs
module load conda
conda activate wgd
export PATH="$HOME"/git/i-ADHoRe/bin/:"$PATH"
export NUMEXPR_MAX_THREADS=4

#### The intra-specific synteny inference
/usr/bin/time wgd syn -f mRNA "$FAMILIES" <(zcat "$GFF") -ks "$KS" --maxsize 50 --ks_range 0 2 -n 4 -o "$OUTDIR"
