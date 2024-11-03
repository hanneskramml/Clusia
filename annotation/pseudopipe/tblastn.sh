#!/bin/bash
#
#SBATCH --job-name=blast
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=slurm.%x.out
#SBATCH --time=0-01:00:00
##SBATCH --license=scratch-highio

module load ncbiblast

PROT_PREFIX=${1:-split}
DB=${2:-../../input/Cmultiflora.genome.fna}

# split query protein fasta
# seqkit split --by-size "$NSEQ" "$QUERY" -O .
# runtime 3-5 hrs for 1000 seq

PART=$(printf "%04d" ${SLURM_ARRAY_TASK_ID})
FILE="$PROT_PREFIX""$PART"  #split0000
QUERY=$(ls "$FILE")
OUT=output/"$FILE".Out

touch stamps/"$FILE".Start
blastall -p tblastn -m 8 -z 3.1e9 -e .1 -d "$DB" -i "$QUERY" -o "$OUT"
touch stamps/"$FILE".Stamp
