#!/bin/bash

#SBATCH --job-name=wgd.dating
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=0-02:00:00
#SBATCH --output=slurm.%x-%j.out

SPECIES=${1:-"Populus_trichocarpa Salix_purpurea Ricinus_communis Manihot_esculenta Hypericum_perforatum Garcinia_oblongifolia Clusia_multiflora"}
TREE=${2:-timetree.independent.nwk}
SYN=${3:-wgd.syn}
OUTDIR=${4:-wgd.dating}

SPECIES=($SPECIES)
FOCAL=${SPECIES[-1]}

# load/set envs
module load conda
conda activate wgd
export NUMEXPR_MAX_THREADS=1

set -o xtrace

#### The phylogenetic dating of WGDs
if [ ! -d "$OUTDIR"/AnchorKs_FindPeak ]; then
   wgd peak --heuristic wgd.ksd/"$FOCAL".paranome.tsv.ks.tsv -ap "$SYN"/iadhore-out/anchorpoints.txt -sm "$SYN"/iadhore-out/segments.txt -le "$SYN"/iadhore-out/list_elements.txt -mp "$SYN"/iadhore-out/multiplicon_pairs.txt -o "$OUTDIR"
fi

for file in "$OUTDIR"/AnchorKs_FindPeak/*"$FOCAL".paranome.tsv.ks.tsv_95%CI_AP_for_dating_weighted_format.tsv; do
   if [ -z $file ]; then continue; fi
   folder=$(basename $file | sed 's/\..*//')

   if [ ! -f "$OUTDIR"/"$folder"/merge_focus_ap.tsv ]; then wgd dmd -f "$FOCAL" -ap "$file" -n 1 -o "$OUTDIR"/"$folder" -t "$TMPDIR" "${SPECIES[@]}"; fi
   if [ ! -f "$OUTDIR"/"$folder"/mcmctree/Concatenated/pep/mcmc.txt ]; then wgd focus --protcocdating --aamodel lg "$OUTDIR"/"$folder"/merge_focus_ap.tsv -sp "$TREE" -o "$OUTDIR"/"$folder" -d mcmctree -ds 'burnin = 2000' -ds 'sampfreq = 1000' -ds 'nsample = 20000' "${SPECIES[@]}"; fi

   peak=$(echo $folder | sed 's/Peak_\([0-9]\)_.*/\1/')
   awk '{print $(NF-3)}' "$OUTDIR"/"$folder"/mcmctree/Concatenated/pep/mcmc.txt > "$OUTDIR"/"$folder".dates.txt
   python postplot.py postdis "$OUTDIR"/"$folder".dates.txt --percentile 95 --title "WGD peak"$peak" of "$FOCAL --hpd -o "$OUTDIR"/"$folder".WGD_date.pdf
done

#cut -f 8 wgd.dating/Peak2/mcmctree/Concatenated/pep/mcmc.txt > wgd.dating/Peak2.dates.txt
#python postplot.py postdis wgd.dating/Peak2.dates.txt --percentile 95 --title "WGD peak2 of Clusia multiflora" --hpd -o "wgd.dating/Peak2.WGD_date.pdf"
