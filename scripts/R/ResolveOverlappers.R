# This function updates 'unresolved_overlaps' and 'unresolved_genes' fields for revolved genes pair
resolved_pair <-function(gene1, gene2){
  overlappers[gene1, "unresolved_overlaps"] <<- overlappers[gene1, "unresolved_overlaps"] - 1
  overlappers[gene1, "unresolved_genes"] <<- gsub(paste0(", ", gene2, "(, |$)"), "\\1", overlappers[gene1, "unresolved_genes"])
  overlappers[gene1, "unresolved_genes"] <<- gsub(paste0("^", gene2, "(, |$)"), "", overlappers[gene1, "unresolved_genes"])
  overlappers[gene2, "unresolved_overlaps"] <<- overlappers[gene2, "unresolved_overlaps"] - 1
  overlappers[gene2, "unresolved_genes"] <<- gsub(paste0(", ", gene1, "(, |$)"), "\\1", overlappers[gene2, "unresolved_genes"])
  overlappers[gene2, "unresolved_genes"] <<- gsub(paste0("^", gene1, "(, |$)"), "", overlappers[gene2, "unresolved_genes"])
}

resolve_overlaps_pos <- function(gene1, gene2, gene1loc, gene2loc){
  # case gene1_end > gene2_end on positive strand
  if(gene1loc[2] - gene2loc[2] > ends_dist){
    if(overlappers[gene1, "ATrich"] == 1 & overlappers[gene1, "priority"] < overlappers[gene2, "priority"]){
      # Shorten gene2
      overlappers[gene2, "delete_start"] <<- paste(overlappers[gene2, "delete_start"], gene1loc[1])
      overlappers[gene2, "delete_end"] <<- paste(overlappers[gene2, "delete_end"], gene2loc[2])
    }
    else {
      # Shorten gene1
      overlappers[gene1, "delete_start"] <<- paste(overlappers[gene1, "delete_start"], gene1loc[1])
      overlappers[gene1, "delete_end"] <<- paste(overlappers[gene1, "delete_end"], gene2loc[2])
    }
  }
  else{
    if(overlappers[gene1, "priority"] < overlappers[gene2, "priority"]){
      # Shorten gene2
      overlappers[gene2, "delete_start"] <<- paste(overlappers[gene2, "delete_start"], gene1loc[1])
      overlappers[gene2, "delete_end"] <<- paste(overlappers[gene2, "delete_end"], gene2loc[2])
    }
    else if(overlappers[gene1, "priority"] > overlappers[gene2, "priority"]){
      # Shorten gene1
      overlappers[gene1, "delete_start"] <<- paste(overlappers[gene1, "delete_start"], gene1loc[1])
      overlappers[gene1, "delete_end"] <<- paste(overlappers[gene1, "delete_end"], gene2loc[2])
    }
    else {
      # Take a look, pair remains unresolved
    }
  }
  resolved_pair(gene1, gene2)
}

resolve_overlaps_neg <- function(gene1, gene2, gene1loc, gene2loc){
  # case gene1_end > gene2_end on negative strand 
  if(gene1loc[1] - gene2loc[1] > ends_dist){
    if(overlappers[gene2, "ATrich"] == 1 & overlappers[gene1, "priority"] > overlappers[gene2, "priority"]){
      # Shorten gene1
      overlappers[gene1, "delete_start"] <<- paste(overlappers[gene1, "delete_start"], gene1loc[1])
      overlappers[gene1, "delete_end"] <<- paste(overlappers[gene1, "delete_end"], gene2loc[2])
    }
    else {
      # Shorten gene2
      overlappers[gene2, "delete_start"] <<- paste(overlappers[gene2, "delete_start"], gene1loc[1])
      overlappers[gene2, "delete_end"] <<- paste(overlappers[gene2, "delete_end"], gene2loc[2])
    }
  }
  else{
    if(overlappers[gene1, "priority"] > overlappers[gene2, "priority"]){
      # Shorten gene1
      overlappers[gene1, "delete_start"] <<- paste(overlappers[gene1, "delete_start"], gene1loc[1])
      overlappers[gene1, "delete_end"] <<- paste(overlappers[gene1, "delete_end"], gene2loc[2])
    }
    else if(overlappers[gene1, "priority"] < overlappers[gene2, "priority"]){
      # Shorten gene2
      overlappers[gene2, "delete_start"] <<- paste(overlappers[gene2, "delete_start"], gene1loc[1])
      overlappers[gene2, "delete_end"] <<- paste(overlappers[gene2, "delete_end"], gene2loc[2])
    }
    else {
      # Take a look, gene remains unresolved
    }
  }
  resolved_pair(gene1, gene2)
}

library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)

gtf_file <- args[1]
overlappers_file <- args[2]
ends_dist <- as.integer(args[3])
output <- args[4]

# Load overlappers file
overlappers <- read.csv(overlappers_file, header = TRUE, row.names=1)
rownames(overlappers) <- overlappers$gene

# Load gtf file:
genome_annotation <- rtracklayer::import(con = gtf_file, format = "gtf")
genome_annotation <- as.data.frame(genome_annotation)
genome_annotation <- genome_annotation[genome_annotation$gene_id %in% overlappers$gene, ]

genes_df <- genome_annotation[genome_annotation$type == "gene", ]
rownames(genes_df) <- genes_df$gene_id

overlappers$ATrich <- 0
overlappers$delete <- 0
overlappers$delete_start <- ""
overlappers$delete_end <- ""
overlappers$unresolved_overlaps <- overlappers$number_of_gene_overlaps
overlappers$unresolved_genes <- overlappers$overlapping_genes

while(1 %in% overlappers$unresolved_overlaps){
    for(gene1 in overlappers$gene){
        if(overlappers[gene1, "unresolved_overlaps"] != 1) next
        gene2 = overlappers[gene1, "unresolved_genes"]
        strand = genes_df[gene1, "strand"]

        if(strand != genes_df[gene2, "strand"]){
            print(paste(gene1, " and ", gene2, " are on different strands! Skipping this loop."))
            next
        }

        # Gene locations:
        gene1loc = c(genes_df[gene1, "start"], genes_df[gene1, "end"])
        gene2loc = c(genes_df[gene2, "start"], genes_df[gene2, "end"])

        # Gene lengths:
        length1 = gene1loc[2] - gene1loc[1]
        length2 = gene2loc[2] - gene2loc[1]

        # If both genes are low-importance, delete shorter one
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
            # If only one is low-importance, delete that one
        } else if(overlappers[gene1, "priority"] > 20) {
            overlappers[gene1, "delete"] = 1
            resolved_pair(gene1, gene2)
            next
        } else if(overlappers[gene2, "priority"] > 20) {
            overlappers[gene2, "delete"] = 1
            resolved_pair(gene1, gene2)
            next
        }
        # In other cases, take into account priorities and positions (look in the functions for the precise alg)
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

# Save the list with modification for manual inspection (if needed)
write.csv(overlappers, sub("\\.gtf$", ".csv", output))

# Now create updated genome annotation

# Delete those that were marked for deletion
overlappers <- overlappers[overlappers$delete == 0, ]
genome_annotation <- genome_annotation[genome_annotation$gene_id %in% overlappers$gene, ]

# Delete fragments marked for deletion
for(gene in overlappers$gene){
  if(is.na(gene) | gene == "") next
  if(is.na(overlappers[gene, "delete_start"]) | overlappers[gene, "delete_start"] == "") next
  strand = genes_df[gene, "strand"]
  start_dels <- as.integer(strsplit(overlappers[gene, "delete_start"], " ")[[1]])
  start_dels <- start_dels[complete.cases(start_dels)]
  end_dels <- as.integer(strsplit(overlappers[gene, "delete_end"], " ")[[1]])
  end_dels <- end_dels[complete.cases(end_dels)]

  start_del = min(start_dels) - 100
  end_del = max(end_dels) + 100

  gene_entries <- which(genome_annotation$gene_id == gene)

  for(gene_entry in gene_entries){
    if(strand == "+") {
      if((genome_annotation[gene_entry, "start"] >= start_del) & (genome_annotation[gene_entry, "start"] <= end_del)){
      genome_annotation[gene_entry, "start"] = end_del + 1
      }
      else if((genome_annotation[gene_entry, "start"] <= start_del) & (genome_annotation[gene_entry, "end"] >= start_del)){
                genome_annotation[gene_entry, "start"] = end_del + 1
      }
    }
    if(strand == "-") {
      if((genome_annotation[gene_entry, "start"] <= start_del) & (genome_annotation[gene_entry, "end"] >= start_del)){
      genome_annotation[gene_entry, "end"] = start_del - 1
      }
      else if((genome_annotation[gene_entry, "start"] >= start_del) & (genome_annotation[gene_entry, "start"] <= end_del)){
          genome_annotation[gene_entry, "end"] = start_del - 1
      }
    }
  }

 
}
genome_annotation <- genome_annotation[genome_annotation$start <= genome_annotation$end, ]
genome_annotation <- genome_annotation[order(genome_annotation$seqnames, genome_annotation$start), ]

genome_annotation <- GenomicRanges::makeGRangesFromDataFrame(genome_annotation, keep.extra.columns=TRUE, na.rm=TRUE)
rtracklayer::export(genome_annotation, output, format = "gtf")
