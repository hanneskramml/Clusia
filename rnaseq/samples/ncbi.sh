#!/bin/bash

SPECIES=${1:-F}
MAPFILE=${2:-VBCF.sampleIDs.txt}

TAR="Clusia-RNAseq_R13151_20220501.tar.gz"
LANE="HCV5JDSX3_1_R13151_20220501/demultiplexed"

set -o xtrace
mkdir -p rawreads

while read line; do
   VBCFID=$(cut -f1 <<< "$line")
   SID=$(cut -f2 <<< "$line")
   CODE=$(cut -f3 <<< "$line")

   # skip non-matching species
   if [ $(cut -c1 <<< "$CODE") != "$SPECIES" ]; then continue; fi

   echo "*** Processing sample $CODE ***"
   #tar -xzvf "$TAR" "$LANE"/"$VBCFID"/"$VBCFID"_S"$SID"_L001_R?_001.fastq.gz --strip-components 3

   READ1=$(ls "$LANE"/"$VBCFID"/"$VBCFID"_S"$SID"_L001_R1_001.fastq.gz)
   READ2=$(ls "$LANE"/"$VBCFID"/"$VBCFID"_S"$SID"_L001_R2_001.fastq.gz)

   mv "$READ1" rawreads/"$CODE".R1.fastq.gz
   mv "$READ2" rawreads/"$CODE".R2.fastq.gz

done < "${MAPFILE}"
