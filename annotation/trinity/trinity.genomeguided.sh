#!/bin/bash
#
#SBATCH --job-name=trinity.GG
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=slurm.%x.out
#SBATCH --time=4-00:00:00

module load trinityrnaseq

BAM=${1:-Cmultiflora.diploid.alignment.mrna.bam}
OUT=${2:-Cmultiflora.trinity.genomeguided}

# Genome-guided transciptome assembly
Trinity --genome_guided_bam "$BAM" --genome_guided_max_intron 50000 \
      --SS_lib_type RF --CPU 1 --max_memory 50G \
      --grid_exec "../hpc/hpc_cmds_GridRunner.pl --grid_conf ../hpc.conf -c" \
      --workdir "$TMPDIR"/"$OUT" --output "$OUT"
#      --no_distributed_trinity_exec

# /home/apps/trinityrnaseq/2.15.0/util/support_scripts//extract_reads_per_partition.pl --partitions_gff Cmultiflora.diploid.alignment.mrna.bam.norm_200.bam.+.sam.minC1.gff  --coord_sorted_SAM Cmultiflora.diploid.alignment.mrna.bam.norm_200.bam.+.sam --parts_per_directory 2000 --min_reads_per_partition 10  --SS_lib_type RF
# touch Dir_<asdf>.ok
# touch partitions.ok
# sed 's/--workdir/--workdir "$TMPDIR\/trinity"/' trinity_GG.cmds > trinity_GG.hk.cmds
