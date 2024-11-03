#!/bin/bash

conda activate seqkit-2.3.1

python juicebox_assembly_converter.py -a Cmultiflora.phased.scaffolds.k30.juicebox.review2.assembly -f references/Cmultiflora.phased.pcontigs.fna > Cmultiflora.phased.scaffolds.k30.juicebox.review2.agp
grep '^>' Cmultiflora.phased.scaffolds.k30.juicebox.review2.fasta | sed 's/^>\(.*\)/\1/' > Cmultiflora.phased.scaffolds.k30.juicebox.review2.lst

seqkit grep -f Cmultiflora.phased.scaffolds.k30.juicebox.review2.chrom.lst Cmultiflora.phased.scaffolds.k30.juicebox.review2.fasta > Cmultiflora.phased.scaffolds.k30.juicebox.review2.chrom.fna
seqkit grep -f Cmultiflora.phased.scaffolds.k30.juicebox.review2.scaff.lst Cmultiflora.phased.scaffolds.k30.juicebox.review2.fasta > Cmultiflora.phased.scaffolds.k30.juicebox.review2.scaff.fna
seqkit grep -f Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.lst Cmultiflora.phased.scaffolds.k30.juicebox.review2.fasta > Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.fna
seqkit grep -p PGA_scaffold_37__145_contigs__length_2370787 Cmultiflora.phased.scaffolds.k30.juicebox.review2.fasta > Cmultiflora.phased.scaffolds.k30.juicebox.review2.debris.fna

grep -f Cmultiflora.phased.scaffolds.k30.juicebox.review2.scaff.lst Cmultiflora.phased.scaffolds.k30.juicebox.review2.agp > Cmultiflora.phased.scaffolds.k30.juicebox.review2.scaff.agp
grep -f Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.lst Cmultiflora.phased.scaffolds.k30.juicebox.review2.agp > Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.agp
grep PGA_scaffold_37__145_contigs__length_2370787 Cmultiflora.phased.scaffolds.k30.juicebox.review2.agp > Cmultiflora.phased.scaffolds.k30.juicebox.review2.debris.agp

# rename headers
awk 'NR==FNR{ a[$1]=$2; next } /^>/{ id=a[substr($0, 2)]; if (id!=""){ print ">" id; next } } 1' Cmultiflora.phased.scaffolds.k30.juicebox.review2.mapping.tsv Cmultiflora.phased.scaffolds.k30.juicebox.review2.chrom.fna$
awk 'NR==FNR{ a[$1]=$2; next } /^>/{ id=a[substr($0, 2)]; if (id!=""){ print ">" id; next } } 1' Cmultiflora.phased.scaffolds.k30.juicebox.review2.mapping.tsv Cmultiflora.phased.scaffolds.k30.juicebox.review2.scaff.fna$

# extract haplotigs
grep PGA_scaffold_36__19_contigs__length_4166560 Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.agp | grep W | cut -f6 | sed 's/:::/___/' > Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.contigs36.lst
seqkit grep -f Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.contigs36.lst Cmultiflora.phased.scaffolds.k30.juicebox.review2.contigs.fna > Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.contigs36.fna
seqkit grep -f Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.scaff.lst Cmultiflora.phased.scaffolds.k30.juicebox.review2.fasta > Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.scaff.fna
cat Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.scaff.fna Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.contigs36.fna > Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.fna
awk 'NR==FNR{ a[$1]=$2; next } /^>/{ id=a[substr($0, 2)]; if (id!=""){ print ">" id; next } } 1' Cmultiflora.phased.scaffolds.k30.juicebox.review2.mapping.tsv Cmultiflora.phased.scaffolds.k30.juicebox.review2.dupseq.fna > Cmultiflora.dupseq.fna
cat Cmultiflora.dupseq.fna Cmultiflora.phased.haplotigs.fna | seqkit sort --by-length -r - > Cmultiflora_v2.haplotigs.softmasked.fna
seqkit seq -u Cmultiflora_v2.haplotigs.softmasked.fna Cmultiflora_v2.haplotigs.fna

# contig to chromosome mapping
awk 'BEGIN {OFS="\t"} NR==FNR{ a[$1]=$2; next } { id=a[$1]; if (id!="" && $5 == "W"){ print id, $2, $3, $6, $7, $8, $9; next } }' Cmultiflora.phased.scaffolds.k30.juicebox.review2.mapping.tsv Cmultiflora.phased.scaffolds.k30.juicebox.review2.agp | sed 's/:::/___/' > Cmultiflora.phased.scaffolds.k30.juicebox.review2.mapping.chrom_contigs.tsv


samtools faidx Cmultiflora_v2.chromosomes.fna
cut -f1,2 Cmultiflora_v2.chromosomes.fna.fai > Cmultiflora_v2.chromosomes.sizes
