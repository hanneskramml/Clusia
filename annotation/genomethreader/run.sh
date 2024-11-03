#!/bin/bash
#SBATCH --job-name=GenomeThreader-pipeline
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=basic
#SBATCH --constraint=array-1core
#SBATCH --output=slurm.pipeline.%j.out
#SBATCH --time=0-0:20:00

# input parameter
GENOME=${1:-Cmultiflora.genome.fna.gz}
PROTEINS=${2:-proteins.faa.gz}
PREFIX=${3:-Cmultiflora}
NSPLIT=${4:-100}

# load envs
module load conda genometools
conda activate seqkit-2.3.1

#### split fastq files
if [ ! -f "split.done" ]; then
   echo "Splitting protein sequences..."
   seqkit split --by-part "$NSPLIT" "$PROTEINS" -O . #proteins.part_001.faa.gz
   # gt splitfasta might work better

   touch "split.done"
fi

#### submit genomethreader jobs
if [ ! -f "submit.done" ]; then
   echo "Submitting genomethreader job array..."
   sbatch -a 1-"$NSPLIT" gth.sh "$GENOME" ${PROTEINS%.faa*} "$PREFIX".gth

   touch "submit.done"
   squeue --me -o "%A %T %j %R %E" | column -t
fi

#### merge and compute consensus spliced alignments
gt gff3 -sort Cmultiflora.gth.part_*.out | gt merge | gt csa | sed 's/gt csa/gth/' - | gt gff3 -sortlines - > Cmultiflora_v2.genes.proteins.gff3
bgzip Cmultiflora_v2.genes.proteins.gff3
tabix -p gff Cmultiflora_v2.genes.proteins.gff3.gz

