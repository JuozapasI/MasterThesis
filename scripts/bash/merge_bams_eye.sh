#!/bin/bash

dir="/tmp/Mazutislab-out/Juozapas/Thesis/data/datasets/eye/"
suffixes=("00" "01" "02" "03" "04" "05" "06")  # Unique prefixes for barcodes

cd ${dir}

# Process each BAM
for i in {1..6}; do
    bam="solo_output.$i/Aligned.sortedByCoord.out.bam"
    suffix=${suffixes[$i]}
    
    # Adjust cell barcodes by adding the suffix
    samtools view -h "$bam" | \
        awk -v sfx="$suffix" 'BEGIN {OFS="\t"} 
            /^@/ {print; next} 
            {for (i=12; i<=NF; i++) if ($i ~ /^CB:Z:/ && $i != "CB:Z:-") $i= $i "_" sfx; print}' | \
        samtools view -b -o "solo_output.$i/modified.bam"
done

# Merge the modified BAMs
samtools merge "merged.bam" solo_output.*/modified.bam

# Sort
samtools sort -@ 8 -o "merged_sorted.bam" "merged.bam" && rm "merged.bam"

# Index the merged BAM
samtools index "merged_sorted.bam"


