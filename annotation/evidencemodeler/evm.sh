#!/bin/bash

module load conda genometools perl genometools gffread samtools
conda activate annotation
#export PATH=/home/user/kramml/git/EVidenceModeler:"$PATH"

# input parameter
PREFIX=${1:-Cmultiflora}
GENOME="$PREFIX".genome.diploid.fna.gz
PREDICTIONS_BRAKER1="$PREFIX".braker1.gff3		# ABINITIO_PREDICTION, Braker1, GFF3-gene-structure
PREDICTIONS_BRAKER2="$PREFIX".braker2.gff3		# ABINITIO_PREDICTION, Braker2, GFF3-gene-structure
PREDICTIONS_PASA="$PREFIX".pasa.orfs.gff3		# OTHER_PREDICTION, transdecoder, GFF3-gene-structure
TRANSCRIPTS_PASA="$PREFIX".pasa.assemblies.gff3		# TRANSCRIPTS, assembler-clusia, GFF3-alignment
TRANSCRIPTS_STRG="$PREFIX".transcripts.gtf		# TRANSCRIPTS, StringTie, GFF3-alignment
PROTEINS_GTH="$PREFIX".genomethreader.gff3		# PROTEIN, GenomeThreader, GFF3-alignment

# intermediate files
PREDICTIONS="$PREFIX".predictions.gff3
TRANSCRIPTS="$PREFIX".transcript_alignments.gff3
PROTEINS="$PREFIX".protein_alignments.gff3
REPEATS="$PREFIX".repeats.gff3

# final output
EVM="$PREFIX".evm.gff3
CDSt="$PREFIX".evm.faa
TRACK="$PREFIX".evm.track.gff3

### prepare input files
if [ "${GENOME##*.}" == "gz" ]; then
   if [ ! -f "${GENOME%.gz}" ]; then pigz -vdfkp4 "$GENOME"; fi
   GENOME="${GENOME%.gz}"
fi

echo "Preparing predictions..."
if [ ! -f "$PREDICTIONS" ]; then
   # { cat "$PREDICTIONS_BRAKER1" | perl ./utils/misc/augustus_GFF3_to_EVM_GFF3.pl - | sed 's/Augustus/Braker1/' -; \
   { cat "$PREDICTIONS_BRAKER1" | sed 's/Liftoff/Braker1/' -; \
      echo; cat "$PREDICTIONS_BRAKER2" | sed 's/AUGUSTUS/Braker2/' -; \
      echo; cat "$PREDICTIONS_PASA" | sed 's/~~/_/'; } > "$PREDICTIONS"
      # echo; perl transdecoder2evm.pl "$PREDICTIONS_PASA" | sed 's/~~/_/'; } > "$PREDICTIONS"

      perl ./utils/gff3_gene_prediction_file_validator.pl "$PREDICTIONS"; echo
fi
ls -la "$PREDICTIONS"

echo "Preparing alignments..."
if [ ! -f "$TRANSCRIPTS" ]; then
   { cat "$TRANSCRIPTS_PASA"; perl ./utils/misc/cufflinks_gtf_to_alignment_gff3.pl "$TRANSCRIPTS_STRG" | sed 's/Cufflinks/StringTie/' -; } > "$TRANSCRIPTS"
fi
if [ ! -f "$PROTEINS" ]; then
   cat "$PROTEINS_GTH" | gt uniq | perl gth2alignment.pl - gth > "$PROTEINS"
fi
ls -la "$TRANSCRIPTS" "$PROTEINS"

echo "Preparing repeats..."
if [ ! -f "$REPEATS" ]; then
   python mask_genome.py "$GENOME" "$PREFIX".gaps.bed "$PREFIX".repeats.bed "$REPEATS"
fi
ls -la "$REPEATS"


### run evidence modeler as slurm job array
echo "Partinioning..."
if [ ! -f "$PREFIX".partitions.list ]; then
   perl ./utils/partition_EVM_inputs.pl --genome "$GENOME" --gene_predictions "$PREDICTIONS" \
      --protein_alignments "$PROTEINS" --transcript_alignments "$TRANSCRIPTS" \ #--repeats "$REPEATS" \
      --segmentSize 1000000 --overlapSize 20000 --partition_listing "$PREFIX".partitions.list
fi
ls -la "$PREFIX".partitions.list

if [ ! -f "$PREFIX".commands.list ]; then
   perl ./utils/write_EVM_commands.pl --partitions "$PREFIX".partitions.list --genome "$GENOME" --gene_predictions "$PREDICTIONS" \
      --protein_alignments "$PROTEINS" --transcript_alignments "$TRANSCRIPTS" \ #--repeats "$REPEATS" \
      --weights /scratch/ecogenomics/clusia/annotation/evidencemodeler/weights.txt --output_file_name evm.out > "$PREFIX".commands.list
fi
ls -la "$PREFIX".commands.list

echo "Submitting slurm array..."
N=$(cat "$PREFIX".commands.list | wc -l)
if [ ! -f "$PREFIX".submit.done ]; then
   sbatch -a 1-"$N" task.sh "$PREFIX"
   touch "$PREFIX".submit.done
   exit 0
else
   K=$(find . -regex ".*evm.out" | wc -l)
   echo "$K/$N jobs processed"
   sleep 3
fi

echo "Creating final outputs..."
if [ ! -f "$EVM" ]; then
   perl ./utils/recombine_EVM_partial_outputs.pl --partitions "$PREFIX".partitions.list --output_file_name evm.out
   perl ./utils/convert_EVM_outputs_to_GFF3.pl --partitions "$PREFIX".partitions.list --output_file_name evm.out --genome "$GENOME"
   find . -maxdepth 2 -regex ".*/evm.out" -exec cat {} \; > "$PREFIX".evm.out
   find . -regex ".*evm.out.gff3" -exec cat {} \; > "$EVM"
fi

if [ ! -f "$CDSt" ]; then gffread "$EVM" -g "$GENOME" -y "$CDSt"; fi

if [ ! -f "$TRACK".gz ]; then
   gt gff3 -sortlines -retainids -tidy "$EVM" > "$TRACK" 2> "$TRACK".log
   bgzip "$TRACK"
   tabix -p gff "$TRACK".gz
fi

ls -la "$EVM" "$CDSt" "$TRACK"*
