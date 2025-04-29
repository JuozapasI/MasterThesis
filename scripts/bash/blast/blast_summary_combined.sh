#!/bin/bash

cd /tmp/Mazutislab-out/Juozapas/Thesis/data/blast/viral/

out1="../blast_counts.txt"
out2="../blast_top_viruses.txt"

AWK_SCRIPT1="../../../scripts/bash/blast_summary.sh"

for file in * ; do
  if [[ -f "$file" ]]; then
    echo "${file%.*}"
    bash "$AWK_SCRIPT1" "$file"
  fi
done > $out1

for file in * ; do
  if [[ -f "$file" ]]; then
    echo "###" "${file%.*}"
    grep "^>" "$file" | sort | uniq -c | sort -k1,1nr | head 
  fi
done > $out2
