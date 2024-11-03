#!/bin/bash
#
#SBATCH --job-name=blast
#SBATCH --cpus-per-task=8
#SBATCH --mem=1G
#SBATCH --output=slurm.%x.out
#SBATCH --time=2-12:00:00
##SBATCH --license=scratch-highio

# ANNOT-pipeline: 3h, 8CPUs, 1G mem
# Cmu=>Ath:       1h, 5CPUs, 1G mem
# Cmu_PH=>sprot   2T

module load ncbiblastplus

QUERY=${1:-query.faa}
DB=${2:-proteindb.faa.gz}
OUT=${3:-${QUERY%.faa}.blastp.tsv}
MAX_SEQS=${4:-1}
PARSE_IDS=${5:-true}	# parse uniprot ids from ncbi accessions (sp|<<uniprot id>>|xxx_yyy)

rsync --copy-links --progress "$QUERY" "$DB" "$TMPDIR"
pushd "$TMPDIR"

if [ "${QUERY##*.}" == "gz" ]; then
   pigz -vdp8 "$QUERY"
   QUERY="${QUERY%.gz}"
fi

if [ "${DB##*.}" == "gz" ]; then
   pigz -vdp8 "$DB"
   DB="${DB%.gz}"
fi

if [ "$PARSE_IDS" = true ]; then
   makeblastdb -in "$DB" -dbtype prot -parse_seqids
else
   makeblastdb -in "$DB" -dbtype prot
fi

blastp -query "$QUERY" -db "$DB" -evalue 1e-5 -max_target_seqs "$MAX_SEQS" -outfmt 6 -num_threads 8 > "$OUT"

popd
rsync --progress "$TMPDIR"/"$OUT" .
