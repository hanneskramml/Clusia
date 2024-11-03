#!/bin/bash
#
# Kmer-based-Subgenome-Mapping
# README: https://github.com/amsession/Kmer-based-Subgenome-Mapping

module load ncbiblastplus perl jellyfish/2.3.0 R

kmer_size=${1:-13}
low_cutoff=${2:-100}
PREFIX=${3:-Cmultiflora}
pairs=${4:-"$PREFIX".homeolog.pairs}
chr_list=${5:-"$PREFIX".chr.list}
chr_leng=${6:-"$PREFIX".chr.lengths}
GENOME=${7:-"$PREFIX".genome.fna.gz}

# zgrep "^>CMU" Cmultiflora.genome.fna.gz | sed 's/>//' > Cmultiflora.chr.list
# bioawk -c fastx '{print $name"\t"length($seq)}' Ghirsutum.genome.fna.gz > Ghirsutum.chr.lengths

OUT="$PREFIX".filt"$low_cutoff".k"$kmer_size"

if [ "${GENOME##*.}" == "gz" ]; then
   if [ ! -f "${GENOME%.gz}" ]; then pigz -vdfkp1 "$GENOME"; fi
   GENOME="${GENOME%.gz}"
fi

if [ ! -f "$GENOME".nhr ]; then
   makeblastdb -dbtype nucl -in "$GENOME" -parse_seqids
fi

kmer_dir="$PREFIX"."$kmer_size"mers
if [ ! -d "$kmer_dir" ]; then
   #hash_size=$(zcat -f "$GENOME" | wc -c | cut -f1 -d ' ')
   hash_size=$(cat "$chr_leng" | cut -f2 | sort -n | tail -1)

   mkdir -p "$kmer_dir" && cd "$kmer_dir"
   perl ../bin/JellyFish_byChr.pl ../"$chr_list" ../"$GENOME" "$hash_size" "$kmer_size"
   cd ..
fi

if [ ! -f "$PREFIX".joined.k"$kmer_size".dump ]; then
   bin/multi_join.sh "$kmer_dir"/sorted.*k"$kmer_size".dump > "$PREFIX".joined.k"$kmer_size".dump
fi

# typically filter to 2 million k-mers or less at this step (use wc -l to check the k-mer count after filtering)
if [ ! -f "$OUT".dump ]; then
   perl bin/filter_JellyRows.pl "$PREFIX".joined.k"$kmer_size".dump "$low_cutoff" > "$OUT".dump
   wc -l "$OUT".dump
fi

if [ ! -f "$OUT".tbl ]; then
   perl bin/chr_list_to_kmer_header.pl "$chr_list" > "$PREFIX".kmer.header
   cat "$PREFIX".kmer.header "$OUT".dump > "$OUT".tbl
fi

if [ ! -f "$OUT".diff.kmer.tbl ]; then
   Rscript bin/Kmap.2subgenome.identification.R "$OUT".tbl "$pairs" "$chr_leng" "$OUT".norm.tbl "$OUT".diff.kmer "$OUT".diff.kmer.tbl #>&Cmultiflora.error.log
   mv Heatmap.pdf "$OUT".heatmap.pdf
   mv Chr.cluster.pdf "$OUT".cluster.pdf
fi

#Rscript bin/Kmap.2subgenome.ANOVA.R Cmultiflora.filt.k13.tbl Cmultiflora.homeolog.pairs Cmultiflora.tukey.out
