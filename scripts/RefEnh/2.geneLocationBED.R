print("2. Making gene location BED file.")

genome_annotation <- rtracklayer::import(con = gtf_file, format = "gtf")
genome_annotation <- as.data.frame(genome_annotation)

gene_ranges_df <- genome_annotation
gene_ranges_df <- gene_ranges_df[gene_ranges_df$type == "gene",]
gene_ranges_df <- GenomicRanges::makeGRangesFromDataFrame(gene_ranges_df, keep.extra.columns=TRUE)
rtracklayer::export(gene_ranges_df, "gene_ranges.gtf", format = "gtf")

rm(genome_annotation)
rm(gene_ranges_df)

## Add "transcript_id """ column to the gtf file to make it compatible with bedtools format
system('awk \'{
       if ($0 ~ "transcript_id") print $0;
       else print $0" transcript_id \"\";"; }\' gene_ranges.gtf > gene_ranges1.gtf')

system('gtf2bed < gene_ranges1.gtf > gene_ranges.bed')

## The following code in R replaces final column with gene name.
gene_ranges = read.table("gene_ranges.bed", sep = "\t")

if(dim(gene_ranges)[1] > 0){
  for (i in 1:dim(gene_ranges)[1])
  {
    a = gene_ranges[i,10]
    res <- stringr::str_match(a, "gene_name\\s*(.*?)\\s*;")
    b = res[,2]
    gene_ranges[i, 10] = b
  }
}

## Save outcome
write.table(gene_ranges, "gene_ranges.bed", sep="\t",row.names=FALSE, col.names=FALSE, quote = FALSE)
system("sort -k1,1 -k2,2n gene_ranges.bed > gene_ranges_sorted.bed")
print("Gene ranges file (gene_ranges_sorted.bed) has been saved in working directory.")

## Remove intermediate files.
system("rm gene_ranges.gtf")
system("rm gene_ranges1.gtf")
system("rm gene_ranges.bed")

print("Done.")