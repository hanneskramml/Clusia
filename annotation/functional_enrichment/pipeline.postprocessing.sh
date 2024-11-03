#!/bin/bash

module load conda gffread genometools bedtools
conda activate agat-1.0.0

# input parameter
PREFIX=${1:-Cmultiflora}
GFF=${2:-"$PREFIX".cid.gff3}
GENOME=${3:-"$PREFIX".genome.diploid.fna.gz}

ANNOT_VERSION=${4:-v2.2}
SPECIES=${5:-Clusia multiflora}


# output
ANNOT_FULL="$PREFIX"_"$ANNOT_VERSION".annotation.gff3
ANNOT_ISO="$PREFIX"_"$ANNOT_VERSION".annotation.longestIsoform.gff3
ANNOT_TRANS="$PREFIX"_"$ANNOT_VERSION".transcripts.fna
ANNOT_CDS="$PREFIX"_"$ANNOT_VERSION".cds.fna
ANNOT_PROT="$PREFIX"_"$ANNOT_VERSION".proteins.faa
ANNOT_TSV="$PREFIX"_"$ANNOT_VERSION".annotation.tsv
ANNOT_STATS="$PREFIX"_"$ANNOT_VERSION".statistics


# postprocessing
if [ "${GFF##*.}" == "gz" ]; then
   if [ ! -f "${GFF%.gz}" ]; then pigz -vdfkp1 "$GFF"; fi
   GFF="${GFF%.gz}"
fi

if [ "${GENOME##*.}" == "gz" ]; then
   if [ ! -f "${GENOME%.gz}" ]; then pigz -vdfkp1 "$GENOME"; fi
   GENOME="${GENOME%.gz}"
fi

ls -la "$GFF" "$GENOME"

if [ ! -f "$ANNOT_FULL" ]; then
   DESC_PREFIX="${GFF%.gff*}".desc
   ATTR="$DESC_PREFIX".attr.gff3
   BED="${ANNOT_FULL%.gff3}".bed

   #agat_convert_sp_gxf2gxf.pl -gff "$GFF" -o "$SORTED"   # sort attributes

   # create file for gene description at level1; TODO: handle cutoff after comma within desc on level2/mRNA
   agat_sp_keep_longest_isoform.pl -gff "$GFF" -o "$PREFIX".tmp
   cat "$PREFIX".tmp | awk -F '\t' '$3=="mRNA" {print $9}' | sed 's/.*Parent=\([^;]*\).*em_desc=\([^;]*\).*/\1\t\2/' | grep -v None | grep -v ";Parent=" > "$DESC_PREFIX".tsv
   { echo -e "ID\tdescription"; cat "$DESC_PREFIX".tsv; } > "$PREFIX".tmp && mv "$PREFIX".tmp "$DESC_PREFIX".tsv
   agat_sq_add_attributes_from_tsv.pl --gff "$GFF" --tsv "$DESC_PREFIX".tsv -o "$DESC_PREFIX".gff3

   agat_sp_manage_attributes.pl --gff "$DESC_PREFIX".gff3 --type mRNA --att em_BRITE,em_COG_cat/cog,em_EC/ec,em_GOs,em_KEGG_Pathway/kegg_pathway,em_KEGG_Reaction/kegg_reaction,em_KEGG_ko/kegg_ko,em_KEGG_Module/kegg_module,em_KEGG_rclass/kegg_rclass,em_KEGG_TC/kegg_tc,em_OGs/og,em_PFAMs/pfam,em_Preferred_name/alt_name,em_desc/desc,em_max_annot_lvl,Ontology_term/ontology,uniprot_id/uniprotID -o "$ATTR"

   { echo "##gff-version 3"; echo "##annot-version $ANNOT_VERSION"; echo "##species $SPECIES"; cat "$ATTR"; } | \
   gt gff3 -sort -tidy -retainids -checkids -addids 2> "$ANNOT_FULL".log > "$ANNOT_FULL"

   agat_convert_sp_gff2bed.pl --gff "$ANNOT_FULL" -o "$BED"
   bedtools sort -i "$BED" > "$PREFIX".tmp && mv "$PREFIX".tmp "$BED"

   createTrack.sh "$ANNOT_FULL"
fi

if [ ! -f "$ANNOT_ISO" ]; then
   BED="${ANNOT_ISO%.gff3}".bed

   agat_sp_keep_longest_isoform.pl -gff "$ANNOT_FULL" -o "$ANNOT_ISO"
   gt gff3 -sort -tidy -retainids -checkids -addids "$ANNOT_ISO" 2> "$ANNOT_ISO".log > tmp && mv tmp "$ANNOT_ISO"

   agat_convert_sp_gff2bed.pl --gff "$ANNOT_ISO" -o "$BED"
   bedtools sort -i "$BED" > tmp && mv tmp "$BED"
fi

if [ ! -f "$ANNOT_PROT" ]; then
   gffread -E -g "$GENOME" -w "$ANNOT_TRANS" -x "$ANNOT_CDS" -y "$ANNOT_PROT" --attrs Name,product,uniprotID,ec "$ANNOT_FULL"
   gffread -E -g "$GENOME" -w "${ANNOT_TRANS%.fna}".primaryTranscript.fna -x "${ANNOT_CDS%.fna}".primaryTranscript.fna -y "${ANNOT_PROT%.faa}".primaryTranscript.faa --attrs Name,product,uniprotID,ec "$ANNOT_ISO"

   sed -i '/^>/ s/;/; /g;s/Name=/Gene=/;s/ec=/EC=/;s/product=/Product=/;s/uniprotID=/UniprotID=/' "$ANNOT_TRANS" "$ANNOT_CDS" "$ANNOT_PROT" "${ANNOT_TRANS%.fna}".primaryTranscript.fna "${ANNOT_CDS%.fna}".primaryTranscript.fna "${ANNOT_PROT%.faa}".primaryTranscript.faa
fi

if [ ! -f "$ANNOT_TSV" ]; then
   gffread -E -g "$GENOME" -o "$ANNOT_TSV" \
      --table @geneid,@id,@chr,@start,@end,@strand,@numexons,@cdslen,Name,alt_name,product,desc,uniprotID,ec,ontology,kegg_ko,pfam,cog,em_target,em_score,em_tcov,em_evalue \
      "$ANNOT_FULL"

   echo -e "Gene\tTranscript\tChr\tStart\tEnd\tStrand\tExonNum\tCDSlen\tName\tAltName\tProduct\tDescription\tUniprotID\tEC\tOntology\tKEGG_ko\tPFAM\tCOG\tem_target\tem_score\tem_scov\tem_evalue" | cat - "$ANNOT_TSV" > tmp && mv tmp "$ANNOT_TSV"
fi

if [ ! -d "$ANNOT_STATS" ]; then
   agat_sp_functional_statistics.pl --gff "$ANNOT_FULL" -o "$ANNOT_STATS"
   agat_sp_statistics.pl --gff "$ANNOT_FULL" -g "$GENOME" -p -o "$ANNOT_STATS"/report_isoforms.txt
fi
