#!/bin/bash

cd /tmp/Mazutislab-out/Juozapas/Thesis/data/blast/viral/plot/

for file in *counts.txt; do
  sample=${file%.*}
  awk -v file=$sample 'BEGIN {OFS = "\t"; split(file, a, "_counts"); sample = a[1]} NR==FNR {top[$2]; next}
       ($2 in top) {species_count[$2] += $1; total += $1} !($2 in top) {species_count["Other"] += $1; total += $1} 
       END {for (sp in species_count) print sample, sp, species_count[sp]/total}' \
      top10_viruses.txt "$file" > "$sample"_proportions.tsv
done


