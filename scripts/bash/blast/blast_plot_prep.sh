#!/bin/bash

cd /tmp/Mazutislab-out/Juozapas/Thesis/data/blast/viral/

mkdir -p plot

for file in * ; do
  if [[ -f "$file" ]]; then
    grep -E '^(Query=|>)' "$file" | \
    awk '/^Query=/ {c = 1; next;} {if (c == 1) {print $0; c = 0;}}' | \
    sort | uniq -c | sort -k1,1nr > plot/${file%.*}_counts.txt
  fi
done 
