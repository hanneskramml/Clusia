#!/bin/bash

## Aligning and QCing Phase Genomics Hi-C Data
## https://phasegenomics.github.io/2019/09/19/hic-alignment-and-qc.html

module load bwa samtools

echo "Assembly file: $1"
echo "HIC reads (fwd): $2"
echo "HIC reads (rev): $3"
echo "Num of threads: $4"

echo "Creating BWA index..."
#bwa index -a bwtsw $1

echo "Creating SAM index..."
#samtools faidx $1

echo "Aligning reads to assembly..."
# local alignment (bwa mem), better use semiglobal for allhic? (bwa aln & bwa sampe)
# removing unmapped, secondary, and supplementary alignments (SAMtools -F 2316), marking (PCR) duplicates (optional)
# bwa mem -5SP -t $4 $1 $2 $3 | ./samblaster | samtools view -S -h -b -F 2316 -@ $4 > aligned.bam
bwa mem -5SP -t $4 $1 $2 $3 | samtools view -S -h -b -F 2316 -@ $4 | samtools sort -@ $4 -o $1.alignment.sorted.bam
