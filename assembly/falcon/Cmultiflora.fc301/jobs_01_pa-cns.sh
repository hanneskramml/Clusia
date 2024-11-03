#!/bin/bash
#
#SBATCH --job-name=C.multi_pa_cns
#SBATCH --cpus-per-task=4
#SBATCH --mem=30000
#SBATCH --constraint="array-4core|array-8core"
#SBATCH --partition=basic
#SBATCH --nice=0
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

ID=$(printf "%05d" ${SLURM_ARRAY_TASK_ID})

cd /scratch/mosys/Clusia_spp_assembly/fc_301/0-rawreads/cns-runs/cns_${ID}/
/bin/bash run-*.bash > run-${SLURM_ARRAY_JOB_ID}_${ID}.stdout 2> run-${SLURM_ARRAY_JOB_ID}_${ID}.stderr
