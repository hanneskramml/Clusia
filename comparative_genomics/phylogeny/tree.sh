#!/bin/bash
#SBATCH --job-name=tree
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=slurm.%x-%j.out
#SBATCH --time=0-00:10:00

module load conda trimal iqtree
conda activate compgen

INPUT=${1:-ITS.fasta}
OUT_PREFIX=${2:-ITS}

muscle -in "$INPUT" -out "$OUT_PREFIX".msa.fna
muscle -in "$OUT_PREFIX".msa.fna -out "$OUT_PREFIX".msa.refined.fna -refine

trimal -in "$OUT_PREFIX".msa.refined.fna -automated1 -phylip -out "$OUT_PREFIX".msa.refined.trimmed.phy
iqtree2 -s "$OUT_PREFIX".msa.refined.trimmed.phy -m GTR+G --alrt 1000 -B 1000

sed -e 's/^/s|/' -e $'s/\t/|/' -e 's/$/|g/' accessions.labels.tsv > accessions.labels.commands.sed
sed -f accessions.labels.commands.sed "$OUT_PREFIX".msa.refined.trimmed.phy.treefile > "$OUT_PREFIX".nw

# esearch -db nucleotide -query "$acc" | efetch -format fasta > $acc.fna
# blastn -query MG779423.fna -db Crosea.genome.fna -evalue 1e-5 -num_threads 4 -outfmt 6 -max_hsps 1 > MG779423.blast.hsp.tsv
# cat MG779423.blast.hsp.tsv | awk -F '\t' '{OFS = FS} {if($10-$9>=0) {print($2,$9-1,$10,$1,$3,"+")} else print($2,$10-1,$9,$1,$3,"-")}' > MG779423.bed
# seqkit subseq --bed MG779423.bed Crosea.genome.fna.gz | seqkit rmdup -s > MG779423.cro.fna
