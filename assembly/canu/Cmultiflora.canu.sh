#!/bin/sh

module unload perl
module load canu/2.1.1 perl

canu -s spec.txt -p Cmultiflora -d Cmultiflora.poly.assembly genomeSize=1.6g -pacbio-raw Cmultiflora.subreads.merged.fna.gz
