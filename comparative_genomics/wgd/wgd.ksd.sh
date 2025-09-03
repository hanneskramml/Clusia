#!/bin/bash

#SBATCH --job-name=wgd.ksd
#SBATCH --cpus-per-task=4
#SBATCH --mem=3G
#SBATCH --time=0-08:00:00
#SBATCH --output=slurm.%x-%j.out

# input parameter
FAMILIES=${1:-wgd.dmd/Clusia_multiflora.paranome.tsv}
CDS=${2:-Clusia_multiflora}
OUTDIR=${3:-wgd.ksd}

# load envs
module load conda
conda activate wgd
export NUMEXPR_MAX_THREADS=4

#### The construction of whole paranome KS age distribution (Cmu: cpu=4, mem=2G, time=8h; Cro: cpu=8, mem=3G, time=19h)
/usr/bin/time wgd ksd "$FAMILIES" "$CDS" -n 4 --cds -t "$TMPDIR" -o "$OUTDIR"

#### The orthologous KS distribution (mem=3G, time=1h for 440/3300 gene families)
#wgd ksd wgd.dmd/global_MRBH.cmu.tsv Clusia_multiflora Garcinia_oblongifolia Hypericum_perforatum Vitis_vinifera -n 4 -t "$TMPDIR" -o "$OUTDIR"
#wgd ksd wgd.dmd/global_MRBH.homoeologs.tsv Clusia_multiflora_H1 Clusia_multiflora_H2 Garcinia_oblongifolia Hypericum_perforatum Vitis_vinifera -n 4 -t "$TMPDIR" -o "$OUTDIR"
