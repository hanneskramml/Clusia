#!/bin/bash
#SBATCH --job-name=genespace
#SBATCH --cpus-per-task=2
#SBATCH --mem=25G
#SBATCH --output=slurm.%x-%j.out
#SBATCH --time=0-03:00:00

# Resource consumption
# eudicots: 4(10) CPUs/tasks, 50G mem, 3h
# clusia_eudicots: 2(10) CPUs/tasks, 30G mem, 1h (2h for pre-diploidized levels of polyploidy)
# clusia.vitis:		2(10) CPUs/tasks, 30G mem, 3h
# clusia.outgroup:	2(10) CPUs/tasks, 25G mem, 3h
# clusia.singlecopy:	2(10) CPUs/tasks, 20G mem, 2h

# do not rerun completed datasets since this leads to erroneous pangene calculation (high mem consumtion) !

source .env
Rscript genespace.R "$1"
