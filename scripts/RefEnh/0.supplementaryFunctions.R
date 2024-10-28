print("0. Loading supplementary functions")

resolved_pair <-function(gene1, gene2){
  overlappers[gene1, "unresolved_overlaps"] <<- overlappers[gene1, "unresolved_overlaps"] - 1
  overlappers[gene1, "unresolved_genes"] <<- gsub(paste0(", ", gene2, "(, |$)"), "\\1", overlappers[gene1, "unresolved_genes"])
  overlappers[gene1, "unresolved_genes"] <<- gsub(paste0("^", gene2, "(, |$)"), "", overlappers[gene1, "unresolved_genes"])
  overlappers[gene2, "unresolved_overlaps"] <<- overlappers[gene2, "unresolved_overlaps"] - 1
  overlappers[gene2, "unresolved_genes"] <<- gsub(paste0(", ", gene1, "(, |$)"), "\\1", overlappers[gene2, "unresolved_genes"])
  overlappers[gene2, "unresolved_genes"] <<- gsub(paste0("^", gene1, "(, |$)"), "", overlappers[gene2, "unresolved_genes"])
}

resolve_embedded <- function(gene1, gene2, gene1loc, gene2loc, strand){
  # gene2 embedded in gene1
  if(overlappers[gene1, "priority"] < overlappers[gene2, "priority"]){
    # Delete gene2
    overlappers[gene1, "classification"] <<- "Keep as is"
    overlappers[gene2, "classification"] <<- "Delete"
    overlappers[gene2, "delete"] <<- 1
  }
  else if(overlappers[gene1, "priority"] == overlappers[gene2, "priority"]){
    if(strand == "+"){
      dist = gene1loc[2] - gene2loc[2]
    }
    else {
      dist = gene2loc[1] - gene1loc[1]
    }
    if(dist > ends_dist){
      # Cut a region in gene1
      overlappers[gene1, "classification"] <<- "Cut"
      overlappers[gene1, "delete_start"] <<- paste(overlappers[gene1, "delete_start"], gene2loc[1])
      overlappers[gene1, "delete_end"] <<- paste(overlappers[gene1, "delete_end"], gene2loc[2])
      overlappers[gene2, "classification"] <<- "Keep as is"
    }
    else{
      # Take a look
      overlappers[gene1, "classification"] <<- "Take a look"
      overlappers[gene2, "classification"] <<- "Take a look"
    }
  }
  else{
    # Cut a region in gene1
    overlappers[gene1, "classification"] <<- "Cut"
    overlappers[gene1, "delete_start"] <<- paste(overlappers[gene1, "delete_start"], gene2loc[1])
    overlappers[gene1, "delete_end"] <<- paste(overlappers[gene1, "delete_end"], gene2loc[2])
    overlappers[gene2, "classification"] <<- "Keep as is"
  }
  resolved_pair(gene1, gene2)
}

resolve_overlaps_pos <- function(gene1, gene2, gene1loc, gene2loc){
  # case gene1_end > gene2_end on positive strand
  if((1-gene1frac)*length1 > ends_dist){
    if(overlappers[gene1, "ATrich"] == 1 & overlappers[gene1, "priority"] < overlappers[gene2, "priority"]){
      # Shorten gene2
      overlappers[gene2, "classification"] <<- "Shorten"
      overlappers[gene2, "delete_start"] <<- paste(overlappers[gene2, "delete_start"], gene1loc[1])
      overlappers[gene2, "delete_end"] <<- paste(overlappers[gene2, "delete_end"], gene2loc[2])
      overlappers[gene1, "classification"] <<- "Keep as is"
    }
    else {
      # Shorten gene1
      overlappers[gene1, "classification"] <<- "Shorten"
      overlappers[gene1, "delete_start"] <<- paste(overlappers[gene1, "delete_start"], gene1loc[1])
      overlappers[gene1, "delete_end"] <<- paste(overlappers[gene1, "delete_end"], gene2loc[2])
      overlappers[gene2, "classification"] <<- "Keep as is"
    }
  }
  else{
    if(overlappers[gene1, "priority"] < overlappers[gene2, "priority"]){
      # Shorten gene2
      overlappers[gene2, "classification"] <<- "Shorten"
      overlappers[gene2, "delete_start"] <<- paste(overlappers[gene2, "delete_start"], gene1loc[1])
      overlappers[gene2, "delete_end"] <<- paste(overlappers[gene2, "delete_end"], gene2loc[2])
      overlappers[gene1, "classification"] <<- "Keep as is"
    }
    else if(overlappers[gene1, "priority"] > overlappers[gene2, "priority"]){
      # Shorten gene1
      overlappers[gene1, "classification"] <<- "Shorten"
      overlappers[gene1, "delete_start"] <<- paste(overlappers[gene1, "delete_start"], gene1loc[1])
      overlappers[gene1, "delete_end"] <<- paste(overlappers[gene1, "delete_end"], gene2loc[2])
      overlappers[gene2, "classification"] <<- "Keep as is"
    }
    else {
      # Take a look
      overlappers[gene1, "classification"] <<- "Take a look"
      overlappers[gene2, "classification"] <<- "Take a look"
    }
  }
  resolved_pair(gene1, gene2)
}

resolve_overlaps_neg <- function(gene1, gene2, gene1loc, gene2loc){
  # case gene1_end > gene2_end on negative strand 
  if((1-gene2frac)*length2 > ends_dist){
    if(overlappers[gene2, "ATrich"] == 1 & overlappers[gene1, "priority"] > overlappers[gene2, "priority"]){
      # Shorten gene1
      overlappers[gene1, "classification"] <<- "Shorten"
      overlappers[gene1, "delete_start"] <<- paste(overlappers[gene1, "delete_start"], gene1loc[1])
      overlappers[gene1, "delete_end"] <<- paste(overlappers[gene1, "delete_end"], gene2loc[2])
      overlappers[gene2, "classification"] <<- "Keep as is"
    }
    else {
      # Shorten gene2
      overlappers[gene2, "classification"] <<- "Shorten"
      overlappers[gene2, "delete_start"] <<- paste(overlappers[gene2, "delete_start"], gene1loc[1])
      overlappers[gene2, "delete_end"] <<- paste(overlappers[gene2, "delete_end"], gene2loc[2])
      overlappers[gene1, "classification"] <<- "Keep as is"
    }
  }
  else{
    if(overlappers[gene1, "priority"] > overlappers[gene2, "priority"]){
      # Shorten gene1
      overlappers[gene1, "classification"] <<- "Shorten"
      overlappers[gene1, "delete_start"] <<- paste(overlappers[gene1, "delete_start"], gene1loc[1])
      overlappers[gene1, "delete_end"] <<- paste(overlappers[gene1, "delete_end"], gene2loc[2])
      overlappers[gene2, "classification"] <<- "Keep as is"
    }
    else if(overlappers[gene1, "priority"] < overlappers[gene2, "priority"]){
      # Shorten gene2
      overlappers[gene2, "classification"] <<- "Shorten"
      overlappers[gene2, "delete_start"] <<- paste(overlappers[gene2, "delete_start"], gene1loc[1])
      overlappers[gene2, "delete_end"] <<- paste(overlappers[gene2, "delete_end"], gene2loc[2])
      overlappers[gene1, "classification"] <<- "Keep as is"
    }
    else {
      # Take a look
      overlappers[gene1, "classification"] <<- "Take a look"
      overlappers[gene2, "classification"] <<- "Take a look"
    }
  }
  resolved_pair(gene1, gene2)
}

# Function to insert new intergenic genes
new_gene <- function(chr, start, end, strand, src, number) {
inter <- data.frame(seqnames = chr,
                    start = start,
                    end = end,
                    width = (end-start+1),
                    strand = strand,
                    source = src,
                    type = "gene",
                    gene_id = paste0("INTERGENIC",number),
                    gene_type = "unknown",
                    gene_name = paste0("INTERGENIC",number),
                    level = 2)
missing_columns <- setdiff(names(genome_annotation), names(inter))
inter[, missing_columns] <- NA
genome_annotation <- rbind(genome_annotation, inter)

inter <- data.frame(seqnames = chr,
                    start = start,
                    end = end,
                    width = (end-start+1),
                    strand = strand,
                    source = src,
                    type = "exon",
                    gene_id = paste0("INTERGENIC",number),
                    gene_type = "unknown",
                    gene_name = paste0("INTERGENIC",number),
                    level = 2,
                    exon_number = 1)
missing_columns <- setdiff(names(genome_annotation), names(inter))
inter[, missing_columns] <- NA
genome_annotation <- rbind(genome_annotation, inter)
return(genome_annotation) 
}

print("Done.")