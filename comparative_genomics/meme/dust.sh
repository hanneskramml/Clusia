#!/bin/bash

module load meme

for f in $(find . -wholename "./fimo/Cmultiflora/promoters/*.fna" -type f); do dust "$f" > "${f%.fna}.dust_masked.fna"; done
for f in $(find . -wholename "./fimo/Cminor/promoters/*.fna" -type f); do dust "$f" > "${f%.fna}.dust_masked.fna"; done
for f in $(find . -wholename "./fimo/Crosea/promoters/*.fna" -type f); do dust "$f" > "${f%.fna}.dust_masked.fna"; done



