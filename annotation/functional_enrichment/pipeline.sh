#!/bin/bash
#
#SBATCH --job-name=ANNOT-pipeline
#SBATCH --cpus-per-task=1
#SBATCH --mem=45G
#SBATCH --output=slurm.%x.out
#SBATCH --time=02:00:00

module load conda gffread genometools bedtools
conda activate annotation

export EGGNOG_DATA_DIR=./db	# global database not available, download with: download_eggnog_data.py --data_dir db

# input parameter
PREFIX=${1:-Cmultiflora}
INPUT=${2:-"$PREFIX".pasa.gff3}
GENOME=${3:-"$PREFIX".genome.diploid.fna.gz}
PROTDB=${4:-proteindb.faa.gz}

ANNOT_VERSION=${5:-v2.2}
GENEID=${6:-Cmu}
SPECIES=${7:-Clusia multiflora}

FILTER_MODELS=${8:-true}
GENERATE_IDS=${9:-true}

# intermediate files
GFF="$PREFIX".sorted.gff3
CDSt="$PREFIX".cds.faa
IPR="$PREFIX".interpro.tsv
BLAST="$PREFIX".blast.tsv
AGAT="$PREFIX".agat.gff3
EMAPPER_SEARCH="$PREFIX".emapper.seed_orthologs
EMAPPER_DECO="$PREFIX".emapper.decorated.gff
HQ_GENES="$PREFIX".HQ_genes.lst
FILTERED="$PREFIX".filtered.gff3
CID="$PREFIX".cid.gff3

# output
# => file namings defined in pipeline.postprocessing.sh - e.g. ANNOT_FULL="$PREFIX"_"$ANNOT_VERSION".annotation.gff3


### Preprocessing
echo "*** Decompressing input files..."
if [ "${INPUT##*.}" == "gz" ]; then
   if [ ! -f "${INPUT%.gz}" ]; then pigz -vdfkp1 "$INPUT"; fi
   INPUT="${INPUT%.gz}"
fi

if [ "${GENOME##*.}" == "gz" ]; then
   if [ ! -f "${GENOME%.gz}" ]; then pigz -vdfkp1 "$GENOME"; fi
   GENOME="${GENOME%.gz}"
fi

if [ "${PROTDB##*.}" == "gz" ]; then
   if [ ! -f "${PROTDB%.gz}" ]; then pigz -vdfkp1 "$PROTDB"; fi
   PROTDB="${PROTDB%.gz}"
fi

ls -la "$INPUT" "$GENOME" "$PROTDB"

echo "*** Sorting genes..."
if [ ! -f "$GFF" ]; then
   gt gff3 -sort -tidy -checkids -retainids -addids "$INPUT" > "$GFF" 2> "$GFF".log
else
   ls -la "$GFF"
fi

echo "*** Extracting protein sequences..."
if [ ! -f "$CDSt" ]; then
   gffread -y "$CDSt" -E -F -V -g "$GENOME" "$GFF"
else
   ls -la "$CDSt"
fi

### Functional annotation
echo "*** Submitting jobs..."
if [ ! -f "$PREFIX".submit.done ]; then
   jid1=$(sbatch --parsable pipeline.interpro.sh "$CDSt" "$IPR")
   jid2=$(sbatch --parsable pipeline.blastp.sh "$CDSt" "$PROTDB" "$BLAST" 10 false)
   jid3=$(sbatch --parsable pipeline.eggnog.sh "$CDSt" "$PREFIX")
   sbatch --dependency=afterok:"$jid1":"$jid2":"$jid3" pipeline.sh "$PREFIX" "$INPUT" "$GENOME" "$PROTDB" "$ANNOT_VERSION" "$GENEID" "$SPECIES" "$FILTER_MODELS" "$GENERATE_IDS"

   touch "$PREFIX".submit.done
   squeue --me -o "%A %T %j %R %E" | column -t
   exit 0
else
   ls -la "$IPR" "$BLAST" "$EMAPPER_SEARCH"
fi

echo "*** Incorporating interpro/blast results..."
if [ ! -f "$AGAT" ]; then
   conda activate agat-1.0.0
   agat_sp_manage_functional_annotation.pl -f "$GFF" -b "$BLAST" --db "$PROTDB" -i "$IPR" --output "$PREFIX".agat
   mv "$PREFIX".agat/"$GFF" "$AGAT"
   rm "$PREFIX".*.agat.log
else
   ls -la "$AGAT"
fi

echo "*** Annotate/incorporate eggnog results..."
if [ ! -f "$EMAPPER_DECO" ]; then
   # 2-step mode (search + annotation). This loads the whole eggnog.db sqlite3 annotation database into memory (~44 GB)
   conda activate eggnog-mapper-2.1.10
   emapper.py -m no_search --annotate_hits_table "$EMAPPER_SEARCH" --tax_scope 33090 --target_taxa 33090 --report_orthologs --excel --decorate_gff "$AGAT" --cpu 1 --dbmem --output "$PREFIX"
else
   ls -la "$EMAPPER_DECO"
fi

### Filer genes based on functional annotation
# 1) keep all PASA models (having transcript evidence)
# 2) keep valid functions / uniprot hits from non-pasa models
# 3) from the rest, keep only those having protein similarity (pident) > 40%, query and subject coverage > 50% based on viridiplantae orthologs (eggnog/emapper)
echo "*** Filtering gene models ..."

if [ ! "$FILTER_MODELS" = true ]; then
      FILTERED="$EMAPPER_DECO"
      echo "SKIPPED!"

elif [ ! -f "$FILTERED" ]; then
   conda activate annotation

   { \
     cat "$EMAPPER_DECO" | awk -F '\t' '$2=="PASA" && $3=="gene" {print $9}' | sed 's/.*ID=\([^;]*\).*/\1/'; \
     cat "$EMAPPER_DECO" | awk -F '\t' '$2=="EVM" && $3=="mRNA" {print $9}' | grep -v "product=hypothetical protein" | sed 's/.*Parent=\([^;]*\).*/\1/'; \
     cat "$EMAPPER_SEARCH" | grep -v "#" | awk '$9>=40 && $10>=50 && $11>=50 {print $1}' | sed 's/\.t.*//'; \
   } | sort -u > "$HQ_GENES"

   python pipeline.extractGeneModels.py "$EMAPPER_DECO" "$HQ_GENES" "$FILTERED"
   #TODO: agat_sp_filter_incomplete_gene_coding_models.pl
else
   ls -la "$FILTERED"
fi

### postprocessing
echo "*** Generating clusia gene IDs..."

if [ ! "$GENERATE_IDS" = true ]; then
      CID="$FILTERED"
      echo "SKIPPED!"

elif [ ! -f "$CID" ]; then
   SORTED="${FILTERED%.gff3}".sorted.gff3
   gt gff3 -sort -tidy -retainids "$FILTERED" > "$SORTED" 2> "$SORTED".log
   python pipeline.generateGeneIDs.py "$SORTED" "$GENEID" "$CID"

else
   ls -la "$CID"
fi

echo "*** Writing final outputs ..."
./pipeline.postprocessing.sh "$PREFIX" "$CID" "$GENOME" "$ANNOT_VERSION" "$SPECIES"

echo "*** Final Output ***"
ls -la "$PREFIX"_"$ANNOT_VERSION".*
