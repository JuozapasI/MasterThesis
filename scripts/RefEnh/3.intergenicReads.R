print("3. Isolating intergenic reads.")

## Firstly, let's take those reads that were not assigned to any genes (i.e. has "GN:Z:-" tag):
seq_data = GenomicAlignments::readGAlignments(bam_file, index=index_file,
                                              param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE),
                                                                              tag = c("GN", "CB", "UB"), what = "flag", tagFilter = list("GN"="-")))
seq_data = data.frame(seq_data)

# Keep only unique reads with valid barcodes:
seq_data$barcodes = paste(seq_data$CB, seq_data$UB, sep="_")
a = nchar(seq_data$barcodes)==28 # logical vector for selecting reads with non-corrupt barcodes
seq_data = seq_data[a,] # exclude all reads that don't have an intact full cellular and molecular barcodes
seq_data = seq_data[!duplicated(seq_data$barcodes),] # exclude all duplicated intergenic reads

# Save extracted intergenic reads as a separate file
seq_data = GenomicRanges::makeGRangesFromDataFrame(seq_data)
seq_data = as(seq_data, "GAlignments")
seq_data = rtracklayer::asBED(seq_data)
rtracklayer::export.bed(seq_data, con = "intergenic_reads1.bed")

# Sort:
system("sort -k1,1 -k2,2n intergenic_reads1.bed > intergenic_reads_sorted1.bed")
system("rm intergenic_reads1.bed")

# Additionally we check and take only those reads that don't intersect with gene bed file,
# as there might be some unassigned reads that actually come from the gene regions.
system("bedtools intersect -v -a intergenic_reads_sorted1.bed -b gene_ranges_sorted.bed > intergenic_reads.bed")

system("sort -k1,1 -k2,2n intergenic_reads.bed > intergenic_reads_sorted.bed")
system("rm intergenic_reads_sorted1.bed")
system("rm intergenic_reads.bed")

rm(seq_data)

print(Done.)