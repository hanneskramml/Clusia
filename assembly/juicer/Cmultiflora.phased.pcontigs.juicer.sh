#!/bin/bash
module load bwa samtools java/1.8u152

./scripts/juicer.sh -z references/Cmultiflora.phased.pcontigs.fna -d Cmultiflora.phased.pcontigs -s MboI -y restriction_sites/Cmultiflora.phased.pcontigs_MboI.txt -e --assembly -q basic -l basic -C 256000000 -Q 12:00:00 -L 2-00:00:00 -D /scratch/ecogenomics/clusia/juicer -S afterdedup
