#!/bin/bash
#Source: https://www.biostars.org/p/112251/#314840

module load bedtools genometools conda
conda activate annotation

# Inout parameters
PREFIX=${1:-Cmultiflora}
GENOME="$PREFIX".genome.fna.gz
GFF="$PREFIX".annotation.gff3.gz
REPEATS="$PREFIX".repeats.bed.gz
#CHROM_SIZE="$PREFIX".genome.size.txt

# Output files
TRANSCRIPTS="$PREFIX".transcripts.bed
INTRONS="$PREFIX".introns.bed
#REPEATS="$PREFIX".repeats.bed
INTRONIC_REPEATS="$PREFIX".introns.repeats.bed
#INTERGENIC="$PREFIX".intergenic.bed
#EXONS="$PREFIX".exons.bed


# BED file with chrom sizes
# awk 'OFS="\t" {print $1, "0", $2}' chromSizes.txt | sort -k1,1 -k2,2n > chromSizes.bed

#?
#LC_ALL=C
#export LC_ALL=C
#sort -k1,1 -k2,2n "$CHROM_SIZE" > "$CHROM_SIZE".sorted
#CHROM_SIZE="$CHROM_SIZE".sorted

# Sort the GFF file
#zcat "$GFF" | awk 'OFS="\t", $1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > "$PREFIX".annotation.sorted.gff3

#?
#gt gff3 -sortlines -retainids -tidy <(zcat "$GFF") > "$PREFIX".annotation.sorted.gff3
#GFF="$PREFIX".annotation.sorted.gff3

# Get intergenic regions
# bedtools complement -i "$GFF" -g "$CHROM_SIZE" > "$INTERGENIC"

# exons
# awk 'OFS="\t", $1 ~ /^#/ {print $0;next} {if ($3 == "exon") print $1, $4-1, $5}' "$GFF" | tr [:blank:] \\t > "$EXONS"

# introns
# bedtools complement -i <(cat "$EXONS" "$INTERGENIC" | sort -k1,1 -k2,2n) -g "$CHROM_SIZE" > "$INTRONS"



# get gene/transcript lengths
#zcat "$GFF" | awk -F '\t' '{OFS = FS} $3=="gene" {split($9,a,"=|;");print $1,($4-1),$5,a[2],0,$7,($5+1-$4)}' > "$PREFIX".genes.bed
zcat "$GFF" | awk -F '\t' '{OFS = FS} $3=="mRNA" {split($9,a,"=|;");print $1,($4-1),$5,a[2],0,$7,($5+1-$4)}' > "$TRANSCRIPTS"
#awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,($3-$2)}' Cmultiflora.annotation.bed > Cmultiflora.transcripts.length.bed

# get intronic regions
python3 getIntronicRegions.py <(zcat "$GFF") "$INTRONS"
awk '{OFS="\t"}{print $0,($3-$2)}' "$INTRONS" > tmp && mv tmp "$INTRONS"

# get repetitive regions
#python3 getRepetitiveRegions.py <(zcat "$GENOME") "$PREFIX".genome.repeats.bed
#cat "$PREFIX".genome.repeats.bed | grep -v $'\t0\t0\t' > "$REPEATS"

# get repeats overlapping in intronic regions
bedtools intersect -a "$INTRONS" -b <(zcat "$REPEATS") -wao > "$INTRONIC_REPEATS"



# cat ../../assembly/juicer/Cmultiflora.phased.scaffolds.k30.juicebox.review2.mapping.chrom_contigs.tsv | awk -F '\t' '{OFS = FS} {print $1,($2-1),$3,$4,0,$7,($3+1-$2)}' > Cmultiflora.contigs.bed
