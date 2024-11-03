#!/bin/bash

module load bedops

REPEATS=$1 #repeat annotation *.out
BED=${REPEATS%.out}.bed

cat "$REPEATS" | rmsk2bed | cut -f-15 > "$BED"
