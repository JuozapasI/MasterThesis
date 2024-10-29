print("7. Resolving overlapping genes.")

# Load gtf file:
genome_annotation <- rtracklayer::import(con = gtf_file, format = "gtf")
genome_annotation <- as.data.frame(genome_annotation)
genome_annotation <- genome_annotation[genome_annotation$type == "gene", ]

# Load overlappers:
overlappers <- read.csv("overlapping_gene_list.csv", row.names = 1)
rownames(overlappers) <- overlappers$gene
genome_annotation <- genome_annotation[genome_annotation$gene_name %in% overlappers$gene, ]
rownames(genome_annotation) <- genome_annotation$gene_name

overlappers$classification <- ""
overlappers$delete <- 0
overlappers$delete_start <- ""
overlappers$delete_end <- ""
overlappers$unresolved_overlaps <- overlappers$number_of_gene_overlaps
overlappers$unresolved_genes <- overlappers$overlapping_genes

overlappers <- overlappers[order(overlappers$number_of_gene_overlaps), ]

loop = 1
# Resolving overlaps (starting only with single overlaps):
while(1 %in% overlappers$unresolved_overlaps){
  print(paste("loop nr ", loop))
  loop = loop + 1
  for(gene1 in overlappers$gene){
    if(overlappers[gene1, "unresolved_overlaps"] != 1) next
    gene2 = overlappers[gene1, "unresolved_genes"]
    
    strand = genome_annotation[gene1, "strand"]
    if(strand != genome_annotation[gene2, "strand"]){
      print(paste(gene1, " and ", gene2, " are on different strands! Skipping this loop."))
      next
    }
    
    # Gene locations:
    gene1loc = c(genome_annotation[gene1, "start"], genome_annotation[gene1, "end"])
    gene2loc = c(genome_annotation[gene2, "start"], genome_annotation[gene2, "end"])
    
    # Gene lengths:
    length1 = gene1loc[2] - gene1loc[1]
    length2 = gene2loc[2] - gene2loc[1]
    
    intersection_size = min(abs(gene1loc[1] - gene2loc[2]), abs(gene1loc[2] - gene2loc[1]),
                            length1, length2)
    
    # Fraction of intersection size compared to gene lengths:
    gene1frac = intersection_size/length1
    gene2frac = intersection_size/length2
    
    if(overlappers[gene1, "priority"] > 20 & overlappers[gene2, "priority"] > 20){
      if(length1 > length2){
        overlappers[gene2, "delete"] = 1
        resolved_pair(gene1, gene2)
        next
      } else {
        overlappers[gene1, "delete"] = 1
        resolved_pair(gene1, gene2)
        next
      }
    } else if(overlappers[gene1, "priority"] > 20) {
      overlappers[gene1, "delete"] = 1
      resolved_pair(gene1, gene2)
      next
    } else if(overlappers[gene2, "priority"] > 20) {
      overlappers[gene2, "delete"] = 1
      resolved_pair(gene1, gene2)
      next
    }
    # Determining the type of intersection:
    if(gene2frac == 1){
      # gene2 embedded in gene1
      resolve_embedded(gene1, gene2, gene1loc, gene2loc, strand)
      next
    }
    if(gene1frac == 1){
      # gene1 embedded in gene2
      resolve_embedded(gene2, gene1, gene2loc, gene1loc, strand)
      next
    }
    if(gene1loc[2] > gene2loc[2] & strand == "+"){
      resolve_overlaps_pos(gene1, gene2, gene1loc, gene2loc)
    }
    else if(gene1loc[2] < gene2loc[2] & strand == "+"){
      resolve_overlaps_pos(gene2, gene1, gene2loc, gene1loc)
    }
    else if(gene1loc[2] > gene2loc[2] & strand == "-"){
      resolve_overlaps_neg(gene1, gene2, gene1loc, gene2loc)
    }
    else {
      resolve_overlaps_neg(gene2, gene1, gene2loc, gene1loc)
    }
  }
}

# Make the table nice
overlappers <- overlappers[order(overlappers$number_of_gene_overlaps, decreasing = TRUE), ]
rownames(overlappers) <- 1:nrow(overlappers)


overlappers$manual_classification <- ""

# Save the results
write.csv(overlappers, "overlapping_gene_list.csv")

print("Done.")