#!/bin/bash
#SBATCH --job-name=hisat
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --partition=basic
#SBATCH --nice=1000
#SBATCH --output=slurm.%x.out
#SBATCH --time=0-01:00:00

module load conda samtools
conda activate hisat2-2.2.1

INDEX=${1:-Cmultiflora.genome}
SID=${SLURM_ARRAY_TASK_ID}

READ1=$(ls reads/*.S"$SID".R1.trimmed.fq.gz)
READ2=$(ls reads/*.S"$SID".R2.trimmed.fq.gz)

PREFIX=$(basename "$READ1")
PREFIX=${PREFIX%.R1.trimmed.fq.gz}

# Annotation
# stranded (RF): https://chipster.csc.fi/manual/library-type-summary.html
hisat2 -p 4 -x "$INDEX" --downstream-transcriptome-assembly --rna-strandness RF -1 "$READ1" -2 "$READ2" --summary-file "$PREFIX".txt | samtools sort -@ 4 > "$PREFIX".bam

# Transcript expression
# hisat2 -p 4 -x "$INDEX" --rna-strandness RF -1 "$READ1" -2 "$READ2" --rg-id "$PREFIX" --rg PL:ILLUMINA --rg PU:HCV5JDSX3.1 --rg LB:LIB-"$PREFIX" --rg SM:"$PREFIX" --summary-file "$PREFIX".txt | samtools sort -@ 4 > "$PREFIX".bam
