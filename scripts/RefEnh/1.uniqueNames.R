print("1. Making all gene names unique.")

# For some reason, there are some gtf file entries with the same name,
# but different gene ids. To avoid confusion and have unique names,
# for those genes we will set gene_name to gene_id.
genome_annotation <- rtracklayer::import(con = original_gtf_file, format = "gtf")
genome_annotation <- as.data.frame(genome_annotation)
genes <- genome_annotation[genome_annotation$type == "gene", ]
non_unique_names <- genes$gene_name[duplicated(genes$gene_name)]
non_unique_names <- unique(non_unique_names)

# Example of non-unique name:
print("Non unique names:")
print(non_unique_names)

# Replace names:
for(gene in non_unique_names){
  genome_annotation$gene_name[which(genome_annotation$gene_name == gene)] <-
    genome_annotation$gene_id[which(genome_annotation$gene_name == gene)]
}

# Save modified gtf file:
genome_annotation = GenomicRanges::makeGRangesFromDataFrame(genome_annotation, keep.extra.columns=TRUE)
rtracklayer::export(genome_annotation, gtf_file, format = "gtf")

print("Done.")