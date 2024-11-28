library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)

gtf_file <- args[1]
extension_candidates_file <- args[2]
new_gene_candidates_file <- args[3]
output <- args[4]

gtf <- rtracklayer::import(con = gtf_file, format = "gtf")
gtf <- as.data.frame(gtf)
extension_candidates <- read.csv(extension_candidates_file, header = FALSE)
colnames(extension_candidates) < c("chrom", "start", "end", "strand", "gene_id")
new_gene_candidates <- read.csv(new_gene_candidates, header = FALSE)
colnames(extension_candidates) < c("seqnames", "start", "end", "strand")
new_gene_candidates$source <- "bam_reads"
new_gene_candidates$feature <- "gene"
new_gene_candidates$score <- "."
new_gene_candidates$frame <- 0

# Extend genes in extension list:
for(i in 1:nrow(extension_candidates)){
    gene = extension_candidates[i, "gene_id"]
    entries = which(gtf$gene_id == gene)
    gene_entry = which(gtf$type[entries] == "gene")
    gene_start = gtf[gene_entry, "start"]
    gene_end = gtf[gene_entry, "end"]

    if(extension_candidates[i, "start"] < gene_start) {extension_limit = extension_candidates[i, "start"]; extend = "start"}
    else {extension_limit = extension_candidates[i, "end"]; extend = "end"}


    for(j in entries){
        if(extend == "start"){
            if(gtf[j, "start"] == gene_start) gtf[j, "start"] <- extension_limit
        }
        else {
            if(gtf[j, "end"] == gene_end) gtf[j,"end"] <- extension_limit
        }
    }
}


# Ensure the same columns
cols_to_add <- setdiff(names(new_gene_candidates), names(gtf))
for (col in cols_to_add) {
  new_gene_candidates[[col]] <- NA
}

# Add new genes
gtf <- rbind(gtf, new_gene_candidates)

# Save modified gft :
gtf <- GenomicRanges::makeGRangesFromDataFrame(gtf, keep.extra.columns=TRUE, na.rm=TRUE)
rtracklayer::export(gtf, output, format = "gtf")
