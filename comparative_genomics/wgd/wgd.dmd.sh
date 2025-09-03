#!/bin/bash

#SBATCH --job-name=wgd.dmd
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=0-00:10:00
#SBATCH --output=slurm.%x-%j.out

CDS=${1:-Clusia_multiflora}

# load/set envs
module load conda
conda activate wgd
export NUMEXPR_MAX_THREADS=4

#### The delineation of whole paranome
wgd dmd "$CDS" --cds -n 4 -o wgd.dmd -t "$TMPDIR"

#### The delineation of global MRBHs
#wgd dmd Clusia_multiflora Garcinia_oblongifolia Hypericum_perforatum Vitis_vinifera --globalmrbh -n 4 -o wgd.dmd.globalmrbh -t "$TMPDIR"
#wgd dmd Clusia_multiflora_H1 Clusia_multiflora_H2 Garcinia_oblongifolia Hypericum_perforatum Vitis_vinifera --globalmrbh -n 4 -o wgd.dmd.globalmrbh
