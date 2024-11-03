#!/bin/bash
#
#SBATCH --job-name=C.minor_pa_cns
#SBATCH --cpus-per-task=8
#SBATCH --mem=60000
##SBATCH --constraint=array-8core
#SBATCH --partition=basic
#SBATCH --nice=1000
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

ID=$(printf "%05d" ${SLURM_ARRAY_TASK_ID})

cd /scratch/mosys/Clusia_spp_assembly/fc_200/0-rawreads/cns-runs/cns_${ID}/
/bin/bash run-*.bash > run-${SLURM_ARRAY_JOB_ID}_${ID}.stdout 2> run-${SLURM_ARRAY_JOB_ID}_${ID}.stderr
