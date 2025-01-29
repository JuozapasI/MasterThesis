#!/bin/bash

int_dir=$1
interval_file=$2 
output_file=$3

> ${output_file}

for file in ${int_dir}/INT*; do
    segment_name=$(basename "$file")
    add_value=$(awk -v segment="$segment_name" '$4 == segment {print $2}' ${interval_file})
    chrom=$(awk -v segment="$segment_name" '$4 == segment {print $1}' ${interval_file})
    grep -v '^#' "$file" | \
    awk -F '\t' -v segment="$segment_name" -v add_value="$add_value" -v chrom="$chrom" \
    'BEGIN {OFS = FS} \
    {$1 = chrom; $2 = $2 "_" segment; $4 = $4 + add_value; $5 = $5 + add_value; print }' >> "$output_file".tmp
done

awk 'BEGIN {FS=OFS="\t"} { 
    $9 = gensub(/transcript_id "g([0-9\.t]+)";/, "transcript_id \"" $2 "." "\\1\";", "g", $9);
    $9 = gensub(/gene_id "g([0-9\.t]+)";/, "gene_id \"" $2 "." "\\1\";", "g", $9);
    if ($9 ~ /^g/) {sub(/^g/, $2 ".", $9);};
    print
}' "$output_file".tmp > "$output_file"

rm "$output_file".tmp
