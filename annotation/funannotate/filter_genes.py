#!/usr/bin/env python

# setup environment first:
# module load conda perl bedtools
# conda activate annotation

import os, sys
from natsort import natsorted
import funannotate.library as lib


# input parameter
GFF = sys.argv[1]     # Cmultiflora.evm.gff3
Genome = sys.argv[2]  # Cmultiflora.genome.diploid.fna

# output
Repeats = Genome + '.repeats.bed'
AssemblyGaps = Genome + '.assembly-gaps.bed'

PREFIX = os.path.splitext(GFF)[0]
Proteins = PREFIX + '.proteins.fa'
BlastModels = PREFIX + '.transposons.tsv'
GFFcleaned = PREFIX + '.cleaned.gff3'

# config
db = './db'
cpus = 2
evalue = 1e-10
pident = 0.5
min_protlen = 50
log_name = PREFIX + '.filter.log'


# create log file
if os.path.isfile(log_name):
    os.remove(log_name)

# initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
lib.log.info('Input GFF: {}; Genome: {}; Output: {}'.format(GFF, Genome, GFFcleaned))

# check that the genome is soft-masked
lib.log.info('Loading genome assembly and parsing soft-masked repetitive sequences')
if not os.path.isfile(Repeats):
    ContigSizes, GenomeLength, maskedSize, percentMask = lib.checkMasklowMem(Genome, Repeats, AssemblyGaps, cpus, tmpdir='tmp')
    lib.log.info('Genome loaded: {:,} scaffolds; {:,} bp; {:.2%} repeats masked'.format(
                len(ContigSizes), GenomeLength, percentMask))

# get protein fasta files
evmCount = lib.countGFFgenes(GFF)
lib.log.info("Generating protein fasta files from {:,} EVM models".format(evmCount))

# translate GFF3 to proteins
if not os.path.isfile(Proteins):
    Genes = {}
    Genes = lib.gff2dict(GFF, Genome, Genes)
    with open(Proteins, 'w') as evmprots:
        for k, v in natsorted(list(Genes.items())):
            for i, x in enumerate(v['ids']):
                Prot = v['protein'][i]
                evmprots.write('>{:} {:}\n{:}\n'.format(x, k, Prot))

# filter bad models
lib.log.info("now filtering out bad gene models (< %i aa in length, transposable elements, etc)." % min_protlen)

if not os.path.isfile(BlastModels):
    lib.RepeatBlast(Proteins, cpus, evalue, db, '.', BlastModels)

if not os.path.isfile(GFFcleaned):
    lib.RemoveBadModels(Proteins, GFF, min_protlen, Repeats, BlastModels, '.', ['overlap', 'blast'], GFFcleaned, pident)

total = lib.countGFFgenes(GFFcleaned)
lib.log.info('{:,} gene models remaining'.format(total))
