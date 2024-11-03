#!/bin/bash

# conda activate hic
# module load samtools
# PREFIX: bwa aln/sampe alignment (excl. file ext .sam, e.g. Cmultiflora.phased.alignment)

GENOME=$1
PREFIX=$2
K=${3:-30}

# add allhic scripts to path
export PATH=/home/user/kramml/git/ALLHiC/scripts/:/home/user/kramml/git/ALLHiC/bin/:$PATH

# preprocess and filter input
PreprocessSAMs.pl $PREFIX.alignment.sam $GENOME MBOI
filterBAM_forHiC.pl $PREFIX.alignment.REduced.paired_only.bam $PREFIX.alignment.REduced.paired_only.clean.sam
samtools view -bt $GENOME.fai -@4 $PREFIX.alignment.REduced.paired_only.clean.sam > $PREFIX.bam

# partition in homologous groups/chromosomes/etc. (ALLHIC_partition covers two functions)
# (1) Extract/counting restriction sites: calculate an empirical distribution of Hi-C link size based on intra-contig links (allhic extract ".$bam." ".$refSeq." --RE ".$esites)
#     Output files: *.clm, *.distribution.txt, *.counts_GATC.txt, *.pairs.txt
# (2) Partition contigs based on prunning bam file (allhic partition $counts_file $pairs_file ".$K." --minREs ".$minRes), differnet default values for minRE (ALLHiC_partition: 25, allhic extract: 10)
#     Output: *.cluster.txt (will be overwritten), *.counts_GATC.[n]g[k].txt

ALLHiC_partition -r $GENOME -b $PREFIX.bam -e MBOI -k $K
# allhic extract Crosea.alignment.allhic_corrected.sorted.unique.REduced.paired_only.bam Crosea.assembly.allhic_corrected.fna --RE GATC
# allhic partition Crosea.alignment.allhic_corrected.sorted.unique.REduced.paired_only.counts_GATC.txt Crosea.alignment.allhic_corrected.sorted.unique.REduced.paired_only.pairs.txt 60 --minREs 10
