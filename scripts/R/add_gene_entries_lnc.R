library(dplyr)
library(rtracklayer)

# Load the GTF file
gtf_file <- "/tmp/Mazutislab-out/Juozapas/Thesis/data/lnc_original.gtf"
output_file <- "/tmp/Mazutislab-out/Juozapas/Thesis/data/lnc.gtf"
gtf <- rtracklayer::import(gtf_file, format = "gtf")
gtf_df <- as.data.frame(gtf)
gtf_df$gene_type <- "lnc"

# Ensure the data contains a gene_id column
if (!"gene_id" %in% colnames(gtf_df)) {
  stop("The GTF file does not contain a 'gene_id' column.")
}

# Create gene entries by grouping by gene_id
gene_entries <- gtf_df %>%
  filter(type == "exon") %>%                        # Only consider exon entries
  group_by(gene_id) %>%                             # Group by gene_id
  summarise(
    seqnames = dplyr::first(seqnames),       # Take the chromosome
    start = min(start),                             # Start is the minimum exon start
    end = max(end),                                 # End is the maximum exon end
    strand = dplyr::first(strand),           # Take the strand information
    score = ".",                                    # Placeholder score
    source = dplyr::first(source),           # Source (e.g., "ENSEMBL")
    frame = ".",                                    # Placeholder frame
    type = "gene",                                  # Set type to "gene"
    gene_type = "lnc"
  )

# Convert back to GTF-compatible format
gene_gtf <- makeGRangesFromDataFrame(gene_entries, keep.extra.columns = TRUE)
final_gtf <- c(gtf, gene_gtf)  # Combine original GTF with new gene entries

# Export the modified GTF
rtracklayer::export(final_gtf, output_file, format = "gtf")

cat("Gene entries added and saved to:", output_file, "\n")
