print("6. Identifying overlapping genes.")

# Extract 'gene' entries from gtf file:
genes_df <- rtracklayer::import(con = gtf_file, format = "gtf")
genes_df <- as.data.frame(genes_df)
genes_df = genes_df[genes_df$type == "gene",1:13]
row.names(genes_df) = 1:nrow(genes_df)

gene_names = genes_df$gene_name
genes_df = GenomicRanges::makeGRangesFromDataFrame(genes_df, keep.extra.columns=T) # convert into granges object


overlapper = rep(FALSE, length(gene_names))
number_of_overlaps = rep(0, length(gene_names))
overlapping_genes = rep("", length(gene_names))

# Count overlaps
for (i in 1:length(gene_names)){
  a = sum(GenomicRanges::countOverlaps(genes_df, genes_df[i]))
  if (a>1){
    overlapper[i] = TRUE
    number_of_overlaps[i] = a-1
    conflict_genes = gene_names[as.logical(GenomicRanges::countOverlaps(genes_df, genes_df[i]))]
    conflict_genes = setdiff(conflict_genes, gene_names[i])
    overlapping_genes[i] = paste(conflict_genes, collapse = ', ')
  }
}


# Make nice table:
overlappers = as.data.frame(cbind(gene_names, number_of_overlaps, overlapping_genes))[overlapper,]
colnames(overlappers) = c("gene", "number_of_gene_overlaps", "overlapping_genes")
overlappers$number_of_gene_overlaps = as.integer(overlappers$number_of_gene_overlaps)

# Sort the list:
o = order(overlappers$number_of_gene_overlaps, decreasing = TRUE)
overlappers = overlappers[o,]

# Load gtf file:
genome_annotation <- rtracklayer::import(con = gtf_file, format = "gtf")
genome_annotation <- as.data.frame(genome_annotation)
genome_annotation <- genome_annotation[genome_annotation$type == "gene", ]

# Load overlappers:
rownames(overlappers) <- overlappers$gene
genome_annotation <- genome_annotation[genome_annotation$gene_name %in% overlappers$gene, ]
rownames(genome_annotation) <- genome_annotation$gene_name

valuable_types = c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene",
                   "TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene", "protein_coding")

# Set priority scores, lower score = higher priority
overlappers$priority <- 0
for(gene in overlappers$gene){
  # By gene type:
  if(genome_annotation[gene, "gene_type"] %in% valuable_types) {overlappers[gene, "priority"] = 0}
  else if((genome_annotation[gene, "gene_type"] == "lncRNA") |
          (genome_annotation[gene, "gene_type"] == "lincRNA")) {overlappers[gene, "priority"] = 10}
  else {overlappers[gene, "priority"] = 20}
  
  # By level:
  # overlappers[gene, "priority"] = as.integer(overlappers[gene, "priority"]) +
  #  as.integer(genome_annotation[gene, "level"])
}

# Checking if 5' ends are AT-rich:
overlappers$ATrich <- 0
at_regions <- data.frame(chr = "", start = "", end = "",
                         gene = overlappers$gene, score = ".", strand = ".")
rownames(at_regions) = at_regions$gene
for(gene in at_regions$gene) {
  at_regions[gene, "strand"] = as.character(genome_annotation[gene, "strand"])
  at_regions[gene, "chr"] = as.character(genome_annotation[gene, "seqnames"])
  if(at_regions[gene, "strand"] == "+") {
    at_regions[gene, "start"] = genome_annotation[gene, "start"]
    at_regions[gene, "end"] = genome_annotation[gene, "start"] + 100
  } else {
    at_regions[gene, "start"] = genome_annotation[gene, "end"] - 100
    at_regions[gene, "end"] = genome_annotation[gene, "end"]
  }
}
write.table(at_regions, "gene_plus_ends.bed", sep="\t",row.names=FALSE, col.names=FALSE, quote = FALSE)

system(paste("bedtools getfasta -fi ", fasta_file,
             "-nameOnly -tab -bed gene_plus_ends.bed -fo sequences.csv"))

sequences <- read.csv("sequences.csv", sep = '\t', header = FALSE, row.names = 1)

for(gene in rownames(sequences)) {
  if(grepl("AAAAAAAAAAAAAAAAAAAA|TTTTTTTTTTTTTTTTTTTT", sequences[gene, "V2"])) {
    overlappers[gene, "ATrich"] = 1
  }
}

rownames(overlappers) <- 1:nrow(overlappers)
write.csv(overlappers, "overlapping_gene_list.csv")

print("Done.")
