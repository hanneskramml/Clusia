#!/bin/bash

module load perl genometools samtools gffread
export PATH=/home/user/kramml/git/PASApipeline-v2.5.2/:/home/user/kramml/git/PASApipeline-v2.5.2/bin:/home/user/kramml/git/pblat-2.5.1:/home/user/kramml/git/fasta-36.3.8i/bin:$PATH

# input parameter
PREFIX=${1:-Cmultiflora}
GENOME=${2:-"$PREFIX".genome.diploid.fna.gz}

TRINITY=${3:-"$PREFIX".trinity.denovo.fasta}
TRINITY_GG=${4:-"$PREFIX".trinity.genomeguided.fasta}
GTF=${5:-"$PREFIX".transcripts.gtf}

CONFIG_ALN=${6:-alignAssembly.config}
CONFIG_ANNOT=${7:-annotAssembly.config}

# file names
TRANSCRIPTS="$PREFIX".transcripts.fna
ACCS="$PREFIX".TDN.accs

# preprocess input files
cat "$TRINITY" "$TRINITY_GG" > "$TRANSCRIPTS"
if [ "${GENOME##*.}" == "gz" ]; then
   if [ ! -f "${GENOME%.gz}" ]; then pigz -vdfkp4 "$GENOME"; fi
   GENOME="${GENOME%.gz}"
fi

# extract accessions / full-length transcripts
./pasa/misc_utilities/accession_extractor.pl < "$TRINITY" > "$ACCS"
#./pasa/scripts/extract_FL_transdecoder_entries.pl Cmultiflora.transcripts.fna.transdecoder.gff3 > "$PREFIX".FL.accs

# cleaning trinity transcripts
if [ ! -f "$TRANSCRIPTS".clean ]; then seqclean "$TRANSCRIPTS" -c4; fi

# PASA alignment assembly
sbatch pasa1.align.sh "$GENOME" "$TRANSCRIPTS" "$GTF" "$ACCS" "$CONFIG_ALN"
