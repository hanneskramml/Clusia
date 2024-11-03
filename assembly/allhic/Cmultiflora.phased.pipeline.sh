#!/bin/bash

GENOME=Cmultiflora.phased.pcontigs.fna
READ1=Cmultiflora.hicreads.R1.fastq.gz
READ2=Cmultiflora.hicreads.R2.fastq.gz

PREFIX=Cmultiflora.phased
CHR_CODE=CMU
K=${1:-30}

./pipeline.sh "$GENOME" "$READ1" "$READ2" "$PREFIX" "$CHR_CODE" "$K"
