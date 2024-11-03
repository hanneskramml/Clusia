#!/bin/bash

PREFIX=${1:-Cmultiflora}
TSS=${2:-"$PREFIX".tss.bed}

# replace transcript id by gene id
sed -i "s/\.t[0-9]//" $TSS

# remove duplicate lines in tss.bed file
echo "*** Removing duplicated lines in $TSS file"
sort $TSS | uniq > "$PREFIX".tss.uniq.bed
len_tss=$(cat $TSS | wc -l)
len_uniq=$(cat "$PREFIX".tss.uniq.bed | wc -l)
removed=$((len_tss-len_uniq))
echo "*** DONE! removed $removed duplicates"


