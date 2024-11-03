#!/bin/bash

module load bedtools

# input parameter
PREFIX=${1:-Cmultiflora}
GENES=${2:-"$PREFIX".cam_genes.txt}
OUT_PREFIX=${3:-"$PREFIX".promoters}

GENOME="$PREFIX".genome.fna.gz  
GSIZE="$PREFIX".genome.size.txt
GFF="$PREFIX".annotation.gff3.gz
PROMOTER_LENGTH=1500

echo "*** Checking input files..."
if ! ls -la "$GENOME" "$GSIZE" "$GFF" "$GENES"; then exit 1; fi

# decompress files
if [ "${GENOME##*.}" == "gz" ]; then
   if [ ! -f "${GENOME%.gz}" ]; then pigz -vdfkp1 "$GENOME"; fi
   GENOME="${GENOME%.gz}"
fi

echo "*** Converting annotation..."
if [ ! -f "$PREFIX".genes.bed ]; then
   #agat_convert_sp_gff2bed.pl --gff "$GFF" --sub gene -o "$PREFIX".genes.bed
   zcat "$GFF" | awk -F '\t' '{OFS = FS} $3=="gene" {split($9,a,"=|;");print $1,($4-1),$5,a[2],0,$7,($5+1-$4)}' > "$PREFIX".genes.bed
else
   ls -la "$PREFIX".genes.bed
fi

echo "*** Collecting genes of interest..."
grep -w -f "$GENES" "$PREFIX".genes.bed > "$PREFIX".genes.subset.bed
ngenes=$(cat "$GENES" | wc -l)
nsubset=$(cat "$PREFIX".genes.subset.bed | wc -l)

echo "$nsubset/$ngenes genes selected"
if [ "$nsubset" -eq 0 ]; then
   echo "ERROR! Please check gene file: $GENES"
else
   echo "*** Extracting promoter sequences..."
   bedtools flank -i "$PREFIX".genes.subset.bed -g "$GSIZE" -l "$PROMOTER_LENGTH" -r 0 -s > "$OUT_PREFIX".bed 
   bedtools getfasta -name -s -fi "$GENOME" -bed "$OUT_PREFIX".bed -fo "$OUT_PREFIX".fna
fi

echo "*** Cleanup..."
rm -vf "$GENOME" "$GENOME".fai *.agat.log "$PREFIX".genes.bed "$PREFIX".genes.subset.bed

npromoters=$(cat "$OUT_PREFIX".fna | grep "^>" | wc -l)
echo "*** DONE! Found $npromoters sequences in file $OUT ***"


# hardmasked genome
# awk -F '\t' '{OFS = FS} {print $1,$2,$3+1,$4}' Cmultiflora.repeats.bed > Cmultiflora.repeats.corrected.bed
# bedtools maskfasta -fi Cmultiflora_v2.genome.fna -bed Cmultiflora.repeats.corrected.bed -fo Cmultiflora_v2.genome.hardmasked.fna

# genes selected from ClusiaDB
# awk -F '\t' '$6=="Clusia_rosea" {print $8}' Clusia.cam_genes.tsv > Crosea.cam_genes.txt

# split multi-fasta file
# seqkit split --by-id Crosea.promoters.fna

# extract start codons
# agat_sp_add_start_and_stop.pl --gff Cmultiflora.annotation.gff3 --fasta Cmultiflora.genome.fna --out Cmultiflora.annotation.start_stop_codons.gff3
# cat Cmultiflora.annotation.start_stop_codons.gff3 | awk -F '\t' '{OFS = FS} $3=="start_codon" {split($9,a,"=|;");print $1,($4-1),$5,a[4],0,$7}' > Cmultiflora.start_codons.bed




# *** find missing promoter ***

# extract column 4 containing gene ids from bed file with extracted promoter sequence information and save it
# cut -f 4 "$OUT_PREFIX".bed > "$PREFIX".successfully_extracted_promoters.bed

# get gene id missing in bed file
# missing_gene_id=$(grep -vf "$PREFIX".successfully_extracted_promoters.bed "$GENES")

# check if missing gene id is in bed subset -> yes
# grep missing_gene_id "$PREFIX".genes.subset.bed > "$PREFIX".missing_gene_id.bed

# try to extract promoter.bed -> empty
# bedtools flank -i "$PREFIX".missing_gene_id.bed -g "$GSIZE" -l "$PROMOTER_LENGTH" -r 0 -s > "$PREFIX".missing_promoter.bed

# manually change "$PREFIX".missing_gene_id.bed so that positions of the promoter region are specified
# for fwd strand: start gene (smaller number) -----> (start gene - 1500)
# for rev strand: start gene (larger number) -----> (start gene + 1500)

# try to extract promoter.fna from manually created "$PREFIX".missing_promoter.bed
# bedtools getfasta -name -s -fi "$GENOME" -bed "$PREFIX".missing_promoter.bed -fo Cminor.missing_promoter.fna
# OUTPUT for Cminor: Feature (CMI008004F:118332-119832) beyond the length of CMI008004F size (118332 bp).  Skipping.
