#!/bin/bash

conda activate hic
# module unload perl

ASSEMBLY=${1:-Cmultiflora.phased.scaffolds.k30.juicebox.assembly}
MND=${2:-Cmultiflora.phased.pcontigs.mnd}

#sed -e 's/$/ 30 - - 30  - - -/' merged30.txt > merged30.links.txt
#bash ../../3d-dna/visualize/run-assembly-visualizer.sh -q 30 -a -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100 -p true -c ../../Cmultiflora.phased.scaffolds.k30.juicebox.assembly merged30.links.txt

# build HiC map
# -a => metadata from inter_hist.m and inter.txt (need to be placed in working dir)
bash 3d-dna/visualize/run-assembly-visualizer.sh -q 1 -a -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000 -p true "$ASSEMBLY" "$MND"

# detect mismachtes
bash 3d-dna/edit/run-mismatch-detector.sh -p true Cmultiflora.phased.scaffolds.k30.juicebox.hic

# annotate repeats
bash 3d-dna/edit/run-coverage-analyzer.sh Cmultiflora.phased.scaffolds.k30.juicebox.hic
