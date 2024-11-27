print("Identifying overlapping genes.")
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)

gtf_file <- args[1]
gene_list <- args[2]
output <- args[3]

# Load gtf and overlapping gene list:
gene_ids <- read.csv(gene_list, header = FALSE)
gene_ids = gene_ids$V1
genes_df <- rtracklayer::import(con = gtf_file, format = "gtf")
genes_df <- as.data.frame(genes_df)
genes_df = genes_df[genes_df$gene_id %in% gene_ids, ]
genes_df = genes_df[genes_df$type == "gene", ]
row.names(genes_df) = 1:nrow(genes_df)

genes_df_GR = GenomicRanges::makeGRangesFromDataFrame(genes_df, keep.extra.columns=T) # convert into granges object
gene_ids <- genes_df$gene_id # reassign gene_ids to have the same order as in genes_df

overlapper = rep(FALSE, length(gene_ids))
number_of_overlaps = rep(0, length(gene_ids))
overlapping_genes = rep("", length(gene_ids))

# Count overlaps
for (i in 1:length(gene_ids)){
  a = sum(GenomicRanges::countOverlaps(genes_df_GR, genes_df_GR[i]))
  if (a>1){
    overlapper[i] = TRUE
    number_of_overlaps[i] = a-1
    conflict_genes = gene_ids[as.logical(GenomicRanges::countOverlaps(genes_df_GR, genes_df_GR[i]))]
    conflict_genes = setdiff(conflict_genes, gene_ids[i])
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


valuable_types = c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene",
                   "TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene", "protein_coding")


rownames(genes_df) <- genes_df$gene_id

# Set priority scores, lower score = higher priority
overlappers$priority <- 0
for(gene in overlappers$gene){
  # By gene type:
  if(genes_df[gene, "gene_type"] %in% valuable_types) {overlappers[gene, "priority"] = 0}
  else if((genes_df[gene, "gene_type"] == "lncRNA") |
          (genes_df[gene, "gene_type"] == "lincRNA")) {overlappers[gene, "priority"] = 10}
  else {overlappers[gene, "priority"] = 20}
  
  # By level:
  if ("level" %in% colnames(genes_df)) {
    overlappers[gene, "priority"] = as.integer(overlappers[gene, "priority"]) +
    as.integer(genes_df[gene, "level"])
  }
}

rownames(overlappers) <- 1:nrow(overlappers)
write.csv(overlappers, output)
