#!/bin/bash

module load bedtools
module load meme
module load conda
# conda activate seqkit-2.5.0 # check right version to load by command applocate seqkit

PREFIX=${1:-Cmultiflora}
TSS=${2:-"$PREFIX".tss.uniq.bed}
OUT_PREFIX=${3:-./promoters/"$PREFIX".promoters}
GENES=${4:-"$PREFIX".cam_genes.txt}

GSIZE="$PREFIX".genome.size.txt
GENOME="$PREFIX".genome.fna.gz
PROMOTER_LENGTH=1500

echo "*** Checking input files..."
if ! ls -la "$GENOME" "$GSIZE" "$TSS" "$GENES"; then exit 1; fi 

# decompress files
if [ "${GENOME##*.}" == "gz" ]; then
  if [ ! -f "${GENOME%.gz}" ]; then pigz -vdfkp1 "$GENOME"; fi
  GENOME="${GENOME%.gz}"
fi

echo "*** Extracting promoter sequences..."
bedtools flank -i "$TSS" -g "$GSIZE" -l "$PROMOTER_LENGTH" -r 0 -s > "$OUT_PREFIX".bed 
bedtools getfasta -name -s -fi "$GENOME" -bed "$OUT_PREFIX".bed -fo "$OUT_PREFIX".fna

# change fasta headers so that it only contains gene ids to prevent selection of additional sequences (e.g. Cmu08.g106 would select Cmu08.g1060 etc.)
sed "s/::.*//g" "$OUT_PREFIX".fna > "$OUT_PREFIX".changed_header.fna
mv "$OUT_PREFIX".changed_header.fna "$OUT_PREFIX".fna 

# select cam genes
echo "*** Selecting genes of interest..."
seqkit grep -f "$GENES" "$OUT_PREFIX".fna > "$OUT_PREFIX".subset.fna

npromoters=$(cat "$OUT_PREFIX".fna | grep "^>" | wc -l)
selected_promoters=$(cat "$OUT_PREFIX".subset.fna | grep "^>" | wc -l)

echo "*** DONE! Found $npromoters sequences in file $OUT_PREFIX.fna, and $selected_promoters in file "$OUT_PREFIX".subset.fna ***"

# mask low complexity regions
echo "*** Masking simple repeats..."
dust "$OUT_PREFIX".subset.fna > "$OUT_PREFIX".masked.fna

# split fasta file by gene id
echo "*** Spliting fasta file..."
seqkit split --by-id --quiet "$OUT_PREFIX".masked.fna --out-dir ./promoters/$PREFIX
# --id-regexp "\b(C[a-z][a-z]\d+\.g\d+)"

# clean up filenames
for file in "$OUT_PREFIX".masked.part_*; do
    mv "$file" "${file/masked.part_/}"
done

echo "*** Cleanup..."
rm -vf "$GENOME" "$GENOME".fai 




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
