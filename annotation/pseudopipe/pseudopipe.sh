#!/bin/sh

module load ncbiblast
conda activate tools2

pwd=$(pwd)

# seqkit split -i --by-id-prefix "" -O . <dna/>Cmultiflora.genome.fna
# zcat Cmultiflora.annotation.gff3.gz | awk -F '\t' '{OFS = FS} $3=="exon" {split($9,a,"=|;");print $1,($4-1),$5,a[2],0,$7,($5+1-$4)}' > Cmultiflora.exons.bed
# grep -f ../Cmultiflora.proteins.highConf.lst Cmultiflora.exons.bed > exons/Cmultiflora.exons.highConf.bed
# awk '{close(f);f=$1}{print > f"_exLocs"}' Cmultiflora.exons.highConf.bed

# Usage: ppipe [output dir] [masked dna dir] [input dna dir] [input pep dir] [exon mask dir] [PBS?]
./app/bin/pseudopipe.sh "$pwd"/output "$pwd"/input/Cmultiflora.genome.fna "$pwd"/input/dna "$pwd"/input/Cmultiflora.proteins.highConf.faa "$pwd"/input/exons

# sbatch -a 0-n tblastn.sh
# touch output/blast/jobs
