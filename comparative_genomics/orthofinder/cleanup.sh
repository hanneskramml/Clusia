#!/bin/bash

DATASET=${1:-pangenome}
cd "$DATASET" || return

function archive {
   echo "Archiving $1..."
   tar -c "$1"/ | pigz -p4 - > "$1".tar.gz
   if [ -f "$1".tar.gz ]; then rm -r "$1"; fi
}

if [ -d Gene_Trees ]; then archive "Gene_Trees"; fi
if [ -d MultipleSequenceAlignments ]; then archive "MultipleSequenceAlignments"; fi
if [ -d Orthogroup_Sequences ]; then archive "Orthogroup_Sequences"; fi
if [ -d Resolved_Gene_Trees ]; then archive "Resolved_Gene_Trees"; fi
if [ -d Phylogenetically_Misplaced_Genes ]; then archive "Phylogenetically_Misplaced_Genes"; fi
if [ -d Putative_Xenologs ]; then archive "Putative_Xenologs"; fi
if [ -d Orthologues ]; then archive "Orthologues"; fi

echo "Cleaning WorkingDirectory..."
rm -rf WorkingDirectory/Alignments_ids
rm -rf WorkingDirectory/Sequences_ids
rm -rf WorkingDirectory/Trees_ids
rm -rf WorkingDirectory/Distances_SpeciesTree
rm -rf WorkingDirectory/SpeciesTrees_ids
rm -rf WorkingDirectory/pickle
