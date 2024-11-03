#!/bin/bash

module load conda
conda activate clusia

# run FIMO for Cmultiflora, Cminor, and Crosea promoter sequences (individual run per sequence and motif)
# Crosea
Cro_background="./background/Crosea.background.txt"
Cro_annotation="Crosea.annotation.tsv"

for motif in $(find . -wholename "./fimo_motifs/*.meme" -type f); 
do 
	for promoter in $(find . -wholename "./promoters/Crosea/*.fna" -type f);
	do 
		bash fimo_annot_v1.sh Crosea "$promoter" "$motif" "$Cro_background" "$Cro_annotation"
	done
done

# Cminor
Cmi_background="./background/Cminor.background.txt"
Cmi_annotation="Crosea.annotation.tsv"

for motif in $(find . -wholename "./fimo_motifs/*.meme" -type f);
do
	for promoter in $(find . -wholename "./promoters/Cminor/*.fna" -type f);
	do
		bash fimo_annot_v1.sh Cminor "$promoter" "$motif" "$Cmi_background" "$Cmi_annotation"
	done
done

# Cmultiflora
Cmu_background="./background/Cmultiflora.background.txt"
Cmu_annotation="Cmultiflora.annotation.tsv"

for motif in $(find . -wholename "./fimo_motifs/*.meme" -type f); do for promoter in $(find . -wholename "./promoters/Cmultiflora/*.fna" -type f); do bash fimo_annot_v2.sh Cmultiflora "$promoter" "$motif" "$Cmu_background" "$Cmu_annotation"; done; done

# Bonferroni-correct p-values, filter significant motif matches and concatenate fimo results
python correct_filter_concatenate_fimo_results.py -d "./fimo" -p Cmultiflora Cminor Crosea
