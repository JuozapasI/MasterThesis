library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)

original_gtf_file <- args[1]
new_gtf_file <- args[2]
output <- args[3]

# Loading gtf files
original_gtf <- rtracklayer::import(con = original_gtf_file, format = "gtf")
original_gtf <- as.data.frame(original_gtf)
new_gtf <- rtracklayer::import(con = new_gtf_file, format = "gtf")
new_gtf <- as.data.frame(new_gtf)

genes_orig <- original_gtf[original_gtf$type == "gene", ]
rownames(genes_orig) <- genes_orig$gene_id
genes_n <- new_gtf[new_gtf$type == "gene", ]
rownames(genes_n) <- genes_n$gene_id

# Extract gene entries
genes_original <- original_gtf[original_gtf$type == "gene", ]
genes_new <- new_gtf[new_gtf$type == "gene", ]

gene_ids <- genes_new$gene_id
gene_ids_original <- genes_original$gene_id

# Convert to GRanges objects
genes_original <- GenomicRanges::makeGRangesFromDataFrame(genes_original, keep.extra.columns=T)
genes_new <- GenomicRanges::makeGRangesFromDataFrame(genes_new, keep.extra.columns=T)

# Count overlaps:

overlapper = rep(FALSE, length(gene_ids))
number_of_overlaps = rep(0, length(gene_ids))
overlapping_genes = rep("", length(gene_ids))


for (i in 1:length(gene_ids)){
  a = sum(GenomicRanges::countOverlaps(genes_original, genes_new[i]))
  if (a>0){
    overlapper[i] = TRUE
    number_of_overlaps[i] = a
    conflict_genes = gene_ids_original[as.logical(GenomicRanges::countOverlaps(genes_original, genes_new[i]))]
    overlapping_genes[i] = paste(conflict_genes, collapse = ', ')
  }
}

# Make nice table:
overlappers = as.data.frame(cbind(gene_ids, number_of_overlaps, overlapping_genes))[overlapper,]
colnames(overlappers) = c("gene", "number_of_gene_overlaps", "overlapping_genes")
rownames(overlappers) = overlappers$gene
overlappers$number_of_gene_overlaps = as.integer(overlappers$number_of_gene_overlaps)

# Sort the list:
o = order(overlappers$number_of_gene_overlaps, decreasing = TRUE)
overlappers = overlappers[o,]

# Resolve overlaps by shortening new gene entries

for(gene in overlappers$gene){
    overlaps <- unlist(strsplit(overlappers[gene, "overlapping_genes"], ",\\s*"))
    strand <- genes_n[gene, "strand"]
    for(gene2 in overlaps){
        start <- genes_orig[gene2, "start"] - 100
        end <- genes_orig[gene2, "end"] + 100
        if(strand == "+") {
            if((genes_n[gene, "start"] >= start) & (genes_n[gene, "start"] <= end)){
            genes_n[gene, "start"] = end + 1
            }
            else if((genes_n[gene, "start"] <= start) & (genes_n[gene, "end"] >= start)){
                genes_n[gene, "start"] = end + 1
            }
          }
        if(strand == "-") {
            if((genes_n[gene, "start"] <= start) & (genes_n[gene, "end"] >= start)){
            genes_n[gene, "end"] = start - 1
            }
            else if((genes_n[gene, "start"] >= start) & (genes_n[gene, "start"] <= end)){
                genes_n[gene, "end"] = start - 1
            }
        }
    }
}
genes_n <- genes_n[genes_n$start <= genes_n$end, ]
gene_ids <- genes_n$gene_id

# Now modify all new entries:
new_gtf <- new_gtf[new_gtf$gene_id %in% gene_ids, ]
for(i in 1:nrow(new_gtf)){
    new_gtf[i, "start"] = max(new_gtf[i, "start"], genes_n[new_gtf[i, "gene_id"], "start"])
    new_gtf[i, "end"] = min(new_gtf[i, "end"], genes_n[new_gtf[i, "gene_id"], "end"])
}
new_gtf <- new_gtf[new_gtf$start <= new_gtf$end, ]

# Save modified new gft entries:
new_gtf <- GenomicRanges::makeGRangesFromDataFrame(new_gtf, keep.extra.columns=TRUE, na.rm=TRUE)
rtracklayer::export(new_gtf, output, format = "gtf")