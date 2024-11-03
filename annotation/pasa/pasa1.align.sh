#!/bin/bash
#SBATCH --job-name=pasa
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --output=slurm.%x.out
#SBATCH --time=3-00:00:00

module load gmap conda perl samtools R
conda activate minimap2-2.24
export PATH=/home/user/kramml/git/PASApipeline-v2.5.2/:/home/user/kramml/git/PASApipeline-v2.5.2/bin:/home/user/kramml/git/pblat-2.5.1:$PATH

GENOME=${1:-Cmultiflora.genome.diploid.fna}
TRANSCRIPTS=${2:-Cmultiflora.transcripts.fna}
GTF=${3:-Cmultiflora.transcripts.gtf}
ACCS=${4:-Cmultiflora.TDN.accs}
CONFIG=${5:-alignAssembly.config}

# PASA/Trinity assembly - only aligners really make use of multiple processes
# --TRANSDECODER to predict full length cDNA alignments
# --MAX_INTRON_LENGTH|-I  <int>         (max intron length parameter passed to GMAP or BLAT)  (default: 100000) => defaults to 500000!
# --stringent_alignment_overlap <float>  (suggested: 30.0)  overlapping transcripts must have this min % overlap to be clustered. => default 1bp

Launch_PASA_pipeline.pl \
           -c "$CONFIG" -C -R -g "$GENOME" \
           -t "$TRANSCRIPTS".clean -T -u "$TRANSCRIPTS" \
           --TDN "$ACCS" --trans_gtf "$GTF" --TRANSDECODER \
           --ALIGNERS blat,gmap,minimap2 --transcribed_is_aligned_orient \
           --MAX_INTRON_LENGTH 20000 --stringent_alignment_overlap 30.0 \
           --CPU 4

# EST/Alignment
# Launch_PASA_pipeline.pl \
#           -c "$CONFIG" -C -R -g "$PREFIX".genome.fna \
##           -t "$PREFIX".transcripts.fna.clean -T -u "$PREFIX".transcripts.fna \
#           --trans_gtf "$PREFIX".stringtie.gtf
#           -f "$PREFIX".FL.accs \
#           --ALIGNERS blat,gmap,minimap2 \
#           --CPU 8
