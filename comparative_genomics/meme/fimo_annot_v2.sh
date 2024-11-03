#!/bin/bash

module load meme conda
conda activate clusia

# input
PREFIX=$1
PROMOTERS=$2
MOTIFS=$3
MARKOV_BACKGROUND=$4
ANNOTATION=$5

# define name of output file
motif_filename="$(basename $MOTIFS)"
motif_name="${motif_filename%.meme}"
basename="$(basename $PROMOTERS)"
gene_id_base="${basename%.fna}"
gene_id="${gene_id_base/$PREFIX.promoters./}"
out="./fimo/$PREFIX.$gene_id.$motif_name.fimo"

# run meme/fimo - use thresh = 1 to get total number of tests performed by FIMO 
echo "Scanning sequence "$gene_id" for "$motif_name" motif occurences..."
fimo --norc --thresh 1 --bfile "$MARKOV_BACKGROUND" --oc "$out" --verbosity 1 "$MOTIFS" "$PROMOTERS"




