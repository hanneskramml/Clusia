#!/bin/bash

module load conda genometools
conda activate annotation

if [ "$#" -lt 1 ]; then
   printf "Usage: $(basename "$0") GFF-File [TrackOutputName (*.track.gff3)]\n"
   exit 1
fi

GFF=${1:-annotation.gff3}
TRACK=${2:-${GFF%.gff*}.track.gff3}

gt gff3 -sortlines -retainids -tidy "$GFF" > "$TRACK" 2> "$TRACK".log
bgzip "$TRACK"
tabix -p gff "$TRACK".gz

echo "Output:"
ls -la "$TRACK"*
