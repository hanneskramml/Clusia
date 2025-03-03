#!/bin/bash

module load rsem star

SPECIES=${1:-Cmultiflora}
SPECIES_CODE=${2:-F}
MAPFILE=${3:-samples/VBCF.sampleIDs.txt}

GFF="$SPECIES".annotation.gff3.gz
GENOME="$SPECIES".genome.fna.gz

INDEX_DIR="$SPECIES".rsem.ref
INDEX="$INDEX_DIR"/"$SPECIES"

# build index
if [ ! -f "$INDEX_DIR"/genomeParameters.txt ]; then
   mkdir -p "$INDEX_DIR"
   pigz -dfkvp4 "$GFF" "$GENOME"
   #rsem-prepare-reference -p 4 --hisat2-hca --gff3 "${GFF%.gz}" "${GENOME%.gz}" "$INDEX"
   rsem-prepare-reference -p 32 --star --gff3 "${GFF%.gz}" "${GENOME%.gz}" "$INDEX"
   rm "${GENOME%.gz}" "${GFF%.gz}"
fi

# generate tasks/jobs
while read line; do
   VBCFID=$(cut -f1 <<< "$line")
   SID=$(cut -f2 <<< "$line")
   CODE=$(cut -f3 <<< "$line")

   if [ ! "${CODE:0:1}" == "$SPECIES_CODE" ]; then continue; fi

   PREFIX="$CODE".S"$SID"
   if [ -f samples/"$PREFIX".R1.trimmed.fq.gz ] && [ -f samples/"$PREFIX".R2.trimmed.fq.gz ]; then
      if [ ! "$JOBS" ]; then JOBS="$SID"; else JOBS+=",${SID}"; fi
   fi
done < "${MAPFILE}"

# submit slurm array
sbatch -a "$JOBS" rsem.sh "$INDEX"
