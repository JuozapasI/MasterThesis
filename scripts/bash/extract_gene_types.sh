#!/bin/bash

awk -F'\t' 'BEGIN {OFS = FS} $3 == "gene" { 
    split($9, info, ";");
    for (i in info) {
        if (info[i] ~ /gene_id/) {
            gene_id = substr(info[i], index(info[i], "gene_id") + 9, length(info[i]) - index(info[i], "gene_id") - 9);
        }
        if (info[i] ~ /gene_type/) {
            gene_type = substr(info[i], index(info[i], "gene_type") + 11, length(info[i]) - index(info[i], "gene_type") - 11);
        }
    }
    print gene_id, gene_type;
}' $1

