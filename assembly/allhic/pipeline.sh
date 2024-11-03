#!/bin/bash
#SBATCH --job-name=HiC-pipeline
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=basic
#SBATCH --constraint=array-1core
#SBATCH --output=slurm.pipeline.out
#SBATCH --time=1-00:00:00

# input parameter
GENOME=${1:-Cmultiflora.phased.pcontigs.fna}
READS1=${2:-Cmultiflora.hicreads.R1.fastq.gz}
READS2=${3:-Cmultiflora.hicreads.R2.fastq.gz}
PREFIX=${4:-Cmultiflora.phased}
CHR_CODE=${5:-CMU}
K=${6:-30}

# load envs
eval "$(conda shell.bash hook)"
conda activate hic
module load bwa samtools

# add allhic and juicebox scripts to path
export PATH=/home/user/kramml/git/ALLHiC/scripts/:/home/user/kramml/git/ALLHiC/bin/:/home/user/kramml/git/juicebox_scripts/juicebox_scripts/:$PATH

#### alignment
if [ ! -f "$PREFIX".alignment.sam ]; then
   echo "Generating bwa/sam index..."
   if [ ! -f "$GENOME".sa ]; then bwa index -a bwtsw "$GENOME"; fi
   if [ ! -f "$GENOME".fai ]; then samtools faidx "$GENOME"; fi

   echo "Submitting alignment jobs..."
   jid1=$(sbatch --parsable pipeline.bwa.aln.sh "$GENOME" "$READS1" "$PREFIX".alignment.R1.sai)
   jid2=$(sbatch --parsable pipeline.bwa.aln.sh "$GENOME" "$READS2" "$PREFIX".alignment.R2.sai)
   jid3=$(sbatch --parsable --dependency=afterany:${jid1}:${jid2} pipeline.bwa.sampe.sh "$GENOME" "$PREFIX".alignment.R1.sai "$PREFIX".alignment.R2.sai "$READS1" "$READS2" "$PREFIX".alignment.sam)
   sbatch --dependency=afterok:${jid3} pipeline.sh "$GENOME" "$READS1" "$READS2" "$PREFIX" "$CHR_CODE" "$K"
   exit 0
fi

#### preprocess and filter input
if [ ! -f "$PREFIX".bam ]; then
   echo "Prepocessing alignment..."
   PreprocessSAMs.pl "$PREFIX".alignment.sam "$GENOME" MBOI
   echo "Filtering alignment..."
   filterBAM_forHiC.pl "$PREFIX".alignment.REduced.paired_only.bam "$PREFIX".alignment.REduced.paired_only.clean.sam
   echo "Writing final output ($PREFIX.bam)..."
   samtools view -bt "$GENOME".fai -@4 "$PREFIX".alignment.REduced.paired_only.clean.sam > "$PREFIX".bam
fi

#### partition into homologous groups/chromosomes/etc. (ALLHIC_partition covers two functions)
# (1) Extract/counting restriction sites: calculate an empirical distribution of Hi-C link size based on intra-contig links (allhic extract ".$bam." ".$refSeq." --RE ".$esites)
#     Output files: *.clm, *.distribution.txt, *.counts_GATC.txt, *.pairs.txt
# (2) Partition contigs based on prunning bam file (allhic partition $counts_file $pairs_file ".$K." --minREs ".$minRes), differnet default values for minRE (ALLHiC_partition: 25, allhic extract: 10)
#     Output: *.clusters.txt (will be overwritten), *.counts_GATC.[n]g[k].txt
if [ ! -f "$PREFIX".clusters.txt ]; then
   echo "Counting REsites and partitioning contigs..."
   ALLHiC_partition -r "$GENOME" -b "$PREFIX".bam -e MBOI -k "$K"
fi

#### scaffolding
if [ ! -f "$PREFIX".counts_GATC."$K"g"$K".tour ]; then
   echo "Submitting job array for scaffolding..."
   jid=$(sbatch --parsable -a 1-"$K" pipeline.allhic2.scaffold.sh "$PREFIX" "$K")
   sbatch --dependency=afterok:${jid} pipeline.sh "$GENOME" "$READS1" "$READS2" "$PREFIX" "$CHR_CODE" "$K"
   exit 0
fi

#### build genome (script modified to allow group specific build)
if [ ! -f groups.asm.fasta ]; then
   echo "Build genome..."
   ALLHiC_build "$GENOME" "$PREFIX".counts_GATC."$K"
fi

#### postprocessing
if [ ! -f "$PREFIX".scaffolds.k"$K".fna ]; then
   echo "QC and postprocessing..."
   module load genometools seqkit

   # rename chromosomes
   sed "s/^>"$PREFIX".*."$K"g\(.*\)/>"$CHR_CODE"\1/" groups.asm.fasta > "$PREFIX".scaffolds.k"$K".fna
   sed "s/^"$PREFIX".*."$K"g\(.*\)/"$CHR_CODE"\1/" groups.agp | sort -t $'\t' -k1.4n,1 -k4n,4 > "$PREFIX".scaffolds.k"$K".agp

   # extract unanchored scaffolds
   seqkit grep -rp "$CHR_CODE" groups.asm.fasta > "$PREFIX".scaffolds.k"$K".unanchored.fna

   # QC
   statAGP.pl "$PREFIX".scaffolds.k"$K".agp > "$PREFIX".scaffolds.k"$K".agp.stat
   gt seqstat "$PREFIX".scaffolds.k"$K".fna > "$PREFIX".scaffolds.k"$K".fna.stat
   sbatch pipeline.busco.genome.sh "$PREFIX".scaffolds.k"$K".fna "$PREFIX".scaffolds.k"$K".busco

   # plot
   samtools faidx "$PREFIX".scaffolds.k"$K".fna
   cut -f1,2 "$PREFIX".scaffolds.k"$K".fna.fai | grep -v "^"$CHR_CODE"0" | sort -nr -k2 > "$PREFIX".scaffolds.k"$K".chrom.sizes
   sbatch pipeline.allhic3.plot.sh "$PREFIX" "$K"

   # juicebox
   # agp2assembly.py "$PREFIX".scaffolds.k"$K".agp "$PREFIX".scaffolds.k"$K".juicebox.assembly
fi
