print("4. Generating list of gene extension candidates.")

# Get distances from 3' ends:
system(paste("bedtools closest -a intergenic_reads_sorted.bed -b gene_ranges_sorted.bed -g ",
             genome_index, " -s -D a -id -fu > results.txt"))

# Load distances data:
summary_data = read.table("results.txt", sep = "\t")

# Split into unassigned and assigned reads:
#unassigned_reads <- summary_data[summary_data$V22 == ".", ]
summary_data <- summary_data[summary_data$V22 != ".", ]

distances <- summary_data$V23 * -1
distances <- sort(distances)
hist(distances[distances < 100000], main = "Intergenic Reads",
     xlab = "Distance from 3' ends", ylab = "Count of reads")
x = seq(0, 200000, 1000)
y = sapply(x, function(threshold) sum(distances < threshold)) / length(distances)
plot(x,y, type = "l", main = "Intergenic Reads", 
     xlab = "Distance from 3' ends", ylab = "Fraction of total intergenic reads")



# Leave only the reads within threshold distance:
summary_data = summary_data[summary_data$V23>-1*extensionThreshold, ]

# Group intergenic reads by gene and construct extension candidates data frame:
extension_candidates <- split(summary_data$V23, summary_data$V22)
extension_candidates <- sapply(extension_candidates, length) # counting
extension_candidates <- data.frame(gene = names(extension_candidates),
                                   count = extension_candidates) # convert to data frame
extension_candidates$cpm <- extension_candidates$count/seq_depth*1000000 # add CPM column

# Finding distances between genes.
# Make gene extension candidates bed file:
gene_ranges <- read.table("gene_ranges_sorted.bed", sep = '\t', header = FALSE)
gene_ranges <- gene_ranges[gene_ranges$V10 %in% extension_candidates$gene, ]
write.table(gene_ranges, file = "extension_candidates_unsorted.bed", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
system("sort -k1,1 -k2,2n extension_candidates_unsorted.bed > extension_candidates_sorted.bed")
# Find distances from gene extension candidates 3' ends:
system(paste("bedtools closest -a extension_candidates_sorted.bed -b gene_ranges_sorted.bed -g ",
             genome_index, " -io -fd -s -D a -iu > distances.bed"))

# Read the results:
gene_distances <- read.table("distances.bed", sep = "\t", header = FALSE)

# Add distance column to candidates data frame
extension_candidates <- merge(extension_candidates, gene_distances[, c("V10", "V21")],
                              by.x = "gene", by.y = "V10", all.x = TRUE)

# Sort and make it nice:
extension_candidates <- extension_candidates[order(extension_candidates$count, decreasing = TRUE), ]
names(extension_candidates)[names(extension_candidates) == "V21"] <- "distance"
rownames(extension_candidates) <- seq(1, nrow(extension_candidates))

# Save the list:
write.csv(extension_candidates, "extension_candidates.csv")

# Clean up:
system("rm results.txt")
system("rm distances.bed")
system("rm extension_candidates_unsorted.bed")
system("rm extension_candidates_sorted.bed")

rm(summary_data)

print("Done.")