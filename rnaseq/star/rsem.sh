#!/bin/bash
#SBATCH --job-name=rsem
#SBATCH --cpus-per-task=4
#SBATCH --mem=36GB
#SBATCH --nice=1000
#SBATCH --output=slurm.%x.out
#SBATCH --time=0-03:00:00

module load rsem star

INDEX=${1:-Cmultiflora.rsem.ref/Cmultiflora}
SID=${SLURM_ARRAY_TASK_ID}

READ1=$(ls samples/*.S"$SID".R1.trimmed.fq.gz)
READ2=$(ls samples/*.S"$SID".R2.trimmed.fq.gz)

SAMPLE=$(basename "$READ1")
SAMPLE=${SAMPLE%.R1.trimmed.fq.gz}

rsem-calculate-expression --num-threads 4 --strandedness reverse --star --star-gzipped-read-file --paired-end "$READ1" "$READ2" "$INDEX" "$SAMPLE"   # don't use --temporary-folder "$TMPDIR", STAR/rsem cannot handle this due to differnet physical devices
#rsem-calculate-expression --num-threads 16 --alignments --paired-end "$SID".Aligned.toTranscriptome.out.bam $(pwd)/$GENOME $(pwd)/$SID

