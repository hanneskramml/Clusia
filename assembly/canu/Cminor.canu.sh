#!/bin/sh

module load canu/2.2 perl

canu -s spec.hdusage.txt -p Cminor -d Cminor.reads genomeSize=3g -pacbio-raw Cminor.subreads.merged.fna.gz -stopAfter=trimming
