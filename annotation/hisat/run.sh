#!/bin/bash

module load conda
conda activate hisat2-2.2.1

MAPFILE=${1:-reads/VBCF.sampleIDs.txt}
GENOME=${2:-Cmultiflora.genome.diploid.fna.gz}
SPECIES=${3:-F}

INDEX=${GENOME%.fna.gz}

# build index
if [ ! -f "$INDEX".8.ht2 ]; then
   pigz -dfkvp4 "$GENOME"
   hisat2-build -p4 "${GENOME%.gz}" "$INDEX"
   rm "${GENOME%.gz}"
fi

# generate tasks/jobs
while read line; do
   VBCFID=$(cut -f1 <<< "$line")
   SID=$(cut -f2 <<< "$line")
   CODE=$(cut -f3 <<< "$line")

   if [ ! "${CODE:0:1}" == "$SPECIES" ]; then continue; fi

   PREFIX="$CODE".S"$SID"
   if [ -f reads/"$PREFIX".R1.trimmed.fq.gz ] && [ -f reads/"$PREFIX".R2.trimmed.fq.gz ]; then
      if [ ! "$JOBS" ]; then JOBS="$SID"; else JOBS+=",${SID}"; fi
   fi
done < "${MAPFILE}"

# submit slurm array
sbatch -a "$JOBS" hisat.sh "$INDEX"

# merge bam files after completion:
# samtools merge -@ 4 -o Cmultiflora.diploid.alignment.mrna.bam F*.bam
