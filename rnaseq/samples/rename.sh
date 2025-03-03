#!/bin/bash

MAPFILE=${1:-VBCF.sampleIDs.txt}

while read line; do
   VBCFID=$(cut -f1 <<< "$line")
   SID=$(cut -f2 <<< "$line")
   CODE=$(cut -f3 <<< "$line")

   PREFIX="$CODE".S"$SID"

   if [ -f "$CODE".R1.trimmed.fq.gz ] && [ -f "$CODE".R2.trimmed.fq.gz ]; then
      mv "$CODE".R1.trimmed.fq.gz "$PREFIX".R1.trimmed.fq.gz
      mv "$CODE".R2.trimmed.fq.gz "$PREFIX".R2.trimmed.fq.gz
      mv "$CODE".html "$PREFIX".html
      mv "$CODE".json "$PREFIX".json

      ls "$PREFIX".*
   fi

done < "${MAPFILE}"
