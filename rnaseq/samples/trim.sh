#!/bin/bash

SPECIES=${1:-F}
MAPFILE=${2:-VBCF.sampleIDs.txt}
Q=${3:-30}  #Quality score for trimming
L=${4:-80}  #Read length to accept after trimming

TAR="Clusia-RNAseq_R13151_20220501.tar.gz"
LANE="HCV5JDSX3_1_R13151_20220501/demultiplexed/"

module load fastp/0.23.2

while read line; do
   VBCFID=$(cut -f1 <<< "$line")
   SID=$(cut -f2 <<< "$line")
   CODE=$(cut -f3 <<< "$line")

   # skip non-matching species and samples already processed
   if [ $(cut -c1 <<< "$CODE") != "$SPECIES" ]; then continue; fi
   if [ -f "$CODE".S"$SID".R1.trimmed.fq.gz ] && [ -f "$CODE".S"$SID".R2.trimmed.fq.gz ]; then continue; fi

   echo "*** Processing sample $CODE ***"
   #tar -xzvf "$TAR" "$LANE"/"$VBCFID"/"$VBCFID"_S"$SID"_L001_R?_001.fastq.gz

   READ1=$(ls "$VBCFID"_S"$SID"_L001_R1_001.fastq.gz)
   READ2=$(ls "$VBCFID"_S"$SID"_L001_R2_001.fastq.gz)

   fastp --verbose -q "$Q" -l "$L" --detect_adapter_for_pe --report_title "$CODE (SampleID $SID / $VBCFID)" --html "$CODE".S"$SID".html --json "$CODE".S"$SID".json --thread 4 --in1 "$READ1" --in2 "$READ2" --out1 "$CODE".S"$SID".R1.trimmed.fq.gz --out2 "$CODE".S"$SID".R2.trimmed.fq.gz

   if [ -f "$CODE".S"$SID".R1.trimmed.fq.gz ] && [ -f "$CODE".S"$SID".R2.trimmed.fq.gz ]; then
      rm --verbose "$READ1" "$READ2"
   fi

done < "${MAPFILE}"
