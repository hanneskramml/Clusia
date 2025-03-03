#!/bin/bash
#
#SBATCH --job-name=star
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB
#SBATCH --nice=1000
#SBATCH --output=slurm.%x.out
#SBATCH --time=0-02:00:00

module load star

INDEX=${1:-Cmultiflora.rsem.ref}
SID=${SLURM_ARRAY_TASK_ID}

READ1=$(ls samples/*.S"$SID".R1.trimmed.fq.gz)
READ2=$(ls samples/*.S"$SID".R2.trimmed.fq.gz)

SAMPLE=$(basename "$READ1")
SAMPLE=${SAMPLE%.R1.trimmed.fq.gz}

STAR --genomeDir "$INDEX" --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --runThreadN 4 --genomeLoad NoSharedMemory --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --outSAMheaderHD @HD VN:1.4 SO:unsorted --outFileNamePrefix "$SAMPLE" --readFilesCommand zcat --readFilesIn "$READ1" "$READ2"
