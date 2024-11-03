#!/bin/bash
#SBATCH --job-name=evm
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=slurm.%x.out
#SBATCH --nice=10000
#SBATCH --time=0-01:00:00

module load conda perl
conda activate annotation

PREFIX=${1:-Cmultiflora}

task="$(sed "${SLURM_ARRAY_TASK_ID}q;d" $PREFIX.commands.list)"
eval "$task"
