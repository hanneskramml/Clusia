#!/bin/bash
#
#SBATCH --job-name=RM-pipeline
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --mail-type=end
#SBATCH --output=slurm.pipeline.out
#SBATCH --time=2-00:00:00

# https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
# TE database: DFAM species/clade: rosids to incl. malvids/fabids (famdb.py lineage -ad "rosids") => insufficient families => using plant-specific ensembl db: https://github.com/Ensembl/plant-scripts/releases/download/v0.3/nrTEplantsJune2020.fna.bz2
# serially annotate/mask genomes with RepeatMasker: 1) simple repeats, 2) well-curated repeats from ensemble nrTEplants, 3) classified species-specific repeats, 4) Unknown species-specific repeats

# input parameter
GENOME=${1:-Cmultiflora.haplotigs.fna.gz}
TEdb=${2:-nrTEplantsJune2020.fna.gz}
TElib=${3:-Cmultiflora.TElib.fa.gz}
PREFIX=${4:-Cmultiflora.haplotigs}

# setup environment
module load conda genometools repeatmasker perl
conda activate repeatmasker
conda activate --stack seqkit-2.3.1

# decompress files
if [ "${GENOME##*.}" == "gz" ]; then
   if [ ! -f "${GENOME%.gz}" ]; then pigz -vdfkp1 "$GENOME"; fi
   GENOME="${GENOME%.gz}"
fi

if [ "${TEdb##*.}" == "gz" ]; then
   if [ ! -f "${TEdb%.gz}" ]; then pigz -vdfkp1 "$TEdb"; fi
   TEdb="${TEdb%.gz}"
fi

# preprocessing de-novo library
if [ ! -f "$PREFIX".repeat.families.fna ]; then
   zcat "$TElib" | seqkit fx2tab | awk '{ print "Cmu1_"$0 }' | seqkit tab2fx > "$PREFIX".repeat.families.fna
   cat "$PREFIX".repeat.families.fna | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > "$PREFIX".repeat.families.classified.fna
   cat "$PREFIX".repeat.families.fna | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > "$PREFIX".repeat.families.unknown.fna
fi

# repeatmasker cycles/iterations
REPMASK1="$PREFIX".repeat1.simple
REPMASK2="$PREFIX".repeat2.TEdb
REPMASK3="$PREFIX".repeat3.TElib.classified
REPMASK4="$PREFIX".repeat4.TElib.unknown
REPMASK5="$PREFIX".repeats

# submit repeatmasker jobs to slurm
if [ ! -f "$REPMASK4".fna ]; then
   jid=6938043 # fix dependency issue
   if [ ! -f "$REPMASK1".fna ]; then jid=$(sbatch --parsable pipeline.repmask1.simple.sh "$GENOME" "$REPMASK1"); fi
   if [ ! -f "$REPMASK2".fna ]; then jid=$(sbatch --parsable --dependency=afterok:"$jid" pipeline.repmask2.db.sh "$REPMASK1".fna "$TEdb" "$REPMASK2"); fi
   if [ ! -f "$REPMASK3".fna ]; then jid=$(sbatch --parsable --dependency=afterok:"$jid" pipeline.repmask3.lib.sh "$REPMASK2".fna "$PREFIX".repeat.families.classified.fna "$REPMASK3"); fi
   if [ ! -f "$REPMASK4".fna ]; then jid=$(sbatch --parsable --dependency=afterok:"$jid" pipeline.repmask3.lib.sh "$REPMASK3".fna "$PREFIX".repeat.families.unknown.fna "$REPMASK4"); fi

   sbatch --dependency=afterok:"$jid" pipeline.sh "$GENOME" "$TEdb" "$TElib" "$PREFIX"
   squeue --me -o "%A %T %j %R %E" | column -t

   exit 0
fi

# generate final output - ProcessRepeats requires an unzipped genome
if [ ! -f "$REPMASK5".align ]; then
   cat "$REPMASK1".cat.gz "$REPMASK2".cat.gz "$REPMASK3".cat.gz "$REPMASK4".cat.gz > "$REPMASK5".cat.gz
   ProcessRepeats -a -xsmall -species "Malpighiales" -maskSource "$GENOME" "$REPMASK5".cat.gz
   # rename -v .masked .fna "$GENOME"*
   perl -I /lisc/app/repeatmasker/4.1.6-3.12.1-5.38.2/ rmOutToGFF3.pl "$REPMASK5".out | gt gff3 -sort > "$REPMASK5".gff3
   ./rmOutToBED.sh "$REPMASK5".out
   bgzip "$REPMASK5".gff3
   tabix -p gff "$REPMASK5".gff3.gz
fi

# Visualize repeat landscape
if [ ! -f "$REPMASK5".html ]; then
   size=$(seqkit stats -T "$GENOME" | awk 'NR>1 {print $5}')
   perl -I /lisc/app/repeatmasker/4.1.6-3.12.1-5.38.2/ calcDivergenceFromAlign.pl -s "$REPMASK5".align.divsum "$REPMASK5".align
   ./createRepeatLandscape.pl -div "$REPMASK5".align.divsum -t "$PREFIX" -g "$size" > "$REPMASK5".html
   pigz -vp1 "$REPMASK5".align
fi
