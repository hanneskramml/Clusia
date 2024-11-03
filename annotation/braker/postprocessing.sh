#!/bin/bash

GFF=${1:-Cmultiflora_v2.scaffolds/augustus.hints.gff3}
OUT=${2:-Cmultiflora_v2.genes.abinitio.gff3}
{ echo "##gff-version 3"; cat "$GFF"; } | gt gff3 -sort -retainids > "$OUT"
bgzip "$OUT"
tabix -p gff "$OUT".gz
