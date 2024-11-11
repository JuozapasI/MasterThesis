print("8. Generating final GTF.")

# Load gtf file:
genome_annotation <- rtracklayer::import(con = gtf_file, format = "gtf")
genome_annotation <- as.data.frame(genome_annotation)

# Load overlappers file:
overlappers <- read.csv("overlapping_gene_list.csv")

# Delete genes listed for deletion:
genes_deletion <- overlappers[, c("gene", "delete")]
genes_deletion <- genes_deletion[genes_deletion$delete == 1, ]
genome_annotation <- genome_annotation[!genome_annotation$gene_name %in% genes_deletion$gene, ]

# Delete manually classified genes:
genes_deletion <- overlappers[, c("gene", "manual_classification")]
genes_deletion <- genes_deletion[genes_deletion$manual_classification == "Delete", ]
genome_annotation <- genome_annotation[!genome_annotation$gene_name %in% genes_deletion$gene, ]

# Shorten/delete gene fragments:
genes_shorten <- overlappers[(overlappers$delete_start != "") & ((overlappers$manual_classification == "") | is.na(overlappers$manual_classification)),
                             c("gene", "delete_start", "delete_end")]
genes_shorten <- genes_shorten[genes_shorten$gene %in% genome_annotation$gene_name, ]
rownames(genes_shorten) <- genes_shorten$gene

j = 0
print("Resolving overlapping genes")

for(gene in genes_shorten$gene){
  j = j + 1
  if (j %% 100 == 0) {print(paste(j, " genes resolved."))}
  
  if(is.na(gene) | gene == "") next
  start_del <- as.integer(strsplit(genes_shorten[gene, "delete_start"], " ")[[1]])
  start_del <- start_del[complete.cases(start_del)]
  end_del <- as.integer(strsplit(genes_shorten[gene, "delete_end"], " ")[[1]])
  end_del <- end_del[complete.cases(end_del)]
  for(j in 1:length(start_del)){
    additional_rows <- data.frame()
    gene_entries <- which(genome_annotation$gene_name == gene)
    for(gene_entry in gene_entries){
      if((genome_annotation[gene_entry, "start"] <= start_del[j]) & (genome_annotation[gene_entry, "end"] <= end_del[j])){
        genome_annotation[gene_entry, "end"] = start_del[j] - 1
      }
      else if((genome_annotation[gene_entry, "start"] >= start_del[j]) & (genome_annotation[gene_entry, "end"] >= end_del[j])){
        genome_annotation[gene_entry, "start"] = end_del[j] + 1
      }
      else if((genome_annotation[gene_entry, "start"] >= start_del[j]) & (genome_annotation[gene_entry, "end"] <= end_del[j])){
        genome_annotation[gene_entry, "start"] = -1
      }
      else{
        if (genome_annotation[gene_entry, "type"] != "gene") {
          additional_rows <- rbind(additional_rows, genome_annotation[gene_entry, ])
          genome_annotation[gene_entry, "end"] = start_del[j] - 1
        }
      }
    }
    if(nrow(additional_rows) > 0) {
    additional_rows$start = rep(end_del[j] + 1, nrow(additional_rows))
    genome_annotation <- rbind(genome_annotation, additional_rows)
    }
  }
}
genome_annotation <- genome_annotation[genome_annotation$start != -1, ]
genome_annotation <- genome_annotation[genome_annotation$start <= genome_annotation$end, ]

genome_annotation <- genome_annotation[genome_annotation$start <= genome_annotation$end, ]
genome_annotation <- genome_annotation[order(genome_annotation$seqnames, genome_annotation$start), ]

# rtracklayer::export(genome_annotation, "resolved_overlaps.gtf", format = "gtf")

print("Overlapping genes resolved.")

print("Adding new genes")
new_genes <- read.csv("new_gene_candidates.csv", row.names = 1)

for(i in 1:nrow(new_genes)) genome_annotation <- new_gene(new_genes$chr[i], new_genes$start[i], new_genes$end[i], new_genes$strand[i], "auto", i)
print("New genes added.")

# Extend genes in extension candidates list:
print("Extending genes in extension list.")

extension_candidates <- read.csv("extension_candidates.csv", row.names = 1)
rownames(extension_candidates) <- 1:nrow(extension_candidates)
chr <- read.csv(chr_lengths, sep = '\t', header = FALSE)

j = 0
for (i in 1:nrow(extension_candidates)) {
  j = j + 1
  if (j %% 100 == 0) {print(paste(j, " genes extended."))}
  
  gene = extension_candidates$gene[i]
  if(!(gene %in% genome_annotation$gene_name)) next
  if(is.na(gene) | gene == "") next
  chrm = genome_annotation$seqnames[genome_annotation$gene_name == gene][1]
  chr_length = chr$V2[chr$V1 == chrm][1]
  strd = genome_annotation$strand[genome_annotation$gene_name == gene][1]
  
  # Extend to include the cluster, but not longer than by 10kb
  if(strd == '+') {
    dist = min(10000, extension_candidates$end[i] - genome_annotation$end[tail(which(genome_annotation$gene_name == gene & genome_annotation$type == "gene"), 1)])
    } else if(strd == '-') {
    dist = min(10000, genome_annotation$start[tail(which(genome_annotation$gene_name == gene & genome_annotation$type == "gene"), 1)] - extension_candidates$start[i])
    } else {next}
  
  # Extend gene and exon (last or first) entries, doesn't extending outside the chromosome
  if(strd == '+') {
    genome_annotation$end[tail(which(genome_annotation$gene_name == gene & genome_annotation$type == "exon"), 1)] <-
      min(genome_annotation$end[tail(which(genome_annotation$gene_name == gene & genome_annotation$type == "gene"), 1)] +
            dist, chr_length)
    genome_annotation$end[tail(which(genome_annotation$gene_name == gene & genome_annotation$type == "gene"), 1)] <-
      min(genome_annotation$end[tail(which(genome_annotation$gene_name == gene & genome_annotation$type == "gene"), 1)] +
            dist, chr_length)
  }
  if(strd == '-') {
    genome_annotation$start[which(genome_annotation$gene_name == gene & genome_annotation$type == "exon")[1]] <-
      max(genome_annotation$start[which(genome_annotation$gene_name == gene & genome_annotation$type == "gene")[1]] -
            dist, 0)
    genome_annotation$start[which(genome_annotation$gene_name == gene & genome_annotation$type == "gene")[1]] <-
      max(genome_annotation$start[which(genome_annotation$gene_name == gene & genome_annotation$type == "gene")[1]] -
            dist, 0)
  }
}

print("Extension done.")

genome_annotation <- genome_annotation[genome_annotation$start <= genome_annotation$end, ]
genome_annotation <- genome_annotation[order(genome_annotation$seqnames, genome_annotation$start), ]

genome_annotation <- GenomicRanges::makeGRangesFromDataFrame(genome_annotation, keep.extra.columns=TRUE, na.rm=TRUE)
rtracklayer::export(genome_annotation, "RE.gtf", format = "gtf")

print("Done.")
