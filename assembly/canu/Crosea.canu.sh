#!/bin/sh

module unload perl
module load canu/2.1.1 perl

canu -s spec.txt -p Crosea -d Crosea.poly.assembly genomeSize=3.15g -pacbio-raw Crosea.subreads.merged.fna.gz -gridOptionsBAT="-p himem"
