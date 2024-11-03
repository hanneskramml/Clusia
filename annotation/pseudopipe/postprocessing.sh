#!/bin/bash

module load bedtools

# Generate pseudogene IDs
cat output/pgenes/output_pgenes.txt | sort -k 1,1 -k 2,2n | awk '{printf "%s\tCmu%s.p%s\n", $0, substr($1,4,2), NR-1}' | sed 's/Cmur.p0/id/' > Cmultiflora.pseudogenes.tsv

# Get overlaps with annotated gene models
cat Cmultiflora.pseudogenes.tsv | awk -F '\t' '{OFS = FS} {print $1,$2-1,$3,$15,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,($3+1-$2)}' > Cmultiflora.pseudogenes.bed
bedtools intersect -a <(cat Cmultiflora.pseudogenes.bed | sed '/#chr/d' | awk -F '\t' '{OFS = FS} {print $1,$2,$3,$4,0,$6,($3-$2)}') -b input/Cmultiflora.exons.bed -wo > Cmultiflora.pseudogenes.overlaps.bed
