#!/bin/bash

module load meme

PREFIX=${1:-Cmultiflora}
GENES=${4:-"$PREFIX".non_cam_genes.subset.txt}

# select non cam genes
seqkit grep -r -f "$GENES" ./promoters/"$PREFIX".promoters.fna > ./promoters/"$PREFIX".non_cam.promoters.fna

# create 0-order Markov background model
fasta-get-markov ./promoters/"$PREFIX".non_cam.promoters.fna ./background/"$PREFIX".background.txt