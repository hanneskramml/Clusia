#!/bin/bash
#
#SBATCH --job-name=trinity
#SBATCH --cpus-per-task=6
#SBATCH --mem=50G
#SBATCH --output=slurm.%x.out
#SBATCH --time=8-00:00:00

module load trinityrnaseq

SAMPLES=${1:-Cmultiflora.samples.txt}
OUT=${2:-Cmultiflora.trinity.denovo}

# --prep :Only prepare files (high I/O usage) and stop before kmer counting.

# De novo transcriptome assembly
Trinity --seqType fq --samples_file "$SAMPLES" --SS_lib_type RF --max_chrysalis_cluster_size 50 --CPU 6 --max_memory 50G --grid_exec "../hpc/hpc_cmds_GridRunner.pl --grid_conf ../hpc.conf -c" --workdir "$TMPDIR"/"$OUT" --output "$OUT"
# sed 's/\/tmp\/slurm-5916105/"$TMPDIR"/' recursive_trinity.cmds > recursive_trinity.hk.cmds

# stats
# TrinityStats.pl Cmultiflora.trinity.denovo.fasta > Cmultiflora.trinity.denovo.fasta.stats
