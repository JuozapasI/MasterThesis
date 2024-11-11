print("4. Clustering intergenic reads.")

# Get clusters of unassigned reads:
system("bedtools merge -s -c 6 -o distinct,count < unassigned_reads_sorted.bed > unassigned_reads_clusters.bed")
# Filter to include only large enough clusters:
cpmThresholdCount <- as.integer(seq_depth / 1000000 * cpmThreshold)
system(paste("awk '{if ($5 > ", cpmThresholdCount, ") {print $0}}' unassigned_read_clusters.bed > unassigned_clusters_filtered.bed", sep=""))
# Sort by cluster size:
system("sort -rnk5,5 unassigned_clusters_filtered.bed > unassigned_clusters.bed")

# Sort by coordinate for the bedtools
system("sort -k1,1 -k2,2n unassigned_clusters.bed > unassigned_clusters_sorted.bed")

# Split unassigned clusters into intersecting with genes and non-intersecting (intergenic):
system("bedtools intersect -v -a unassigned_clusters_sorted.bed -b gene_ranges_sorted.bed > intergenic.bed")
system("bedtools intersect -wa -wb -u -a unassigned_clusters_sorted.bed -b gene_ranges_sorted.bed > intersecting.bed")

# Additionally, we check intergenic clusters, if they intersect with the comprehensive annotation:
system("bedtools intersect -v -a intergenic.bed -b gene_ranges_sorted_full.bed > true_intergenic.bed")
system("bedtools intersect -wa -wb -u -a intergenic.bed -b gene_ranges_sorted_full.bed > intergenic_full_annotation.bed")

# For true intergenic 
# Add empty columns to make it compatible with gene_ranges_sorted.bed in bedtools closest:
system("awk 'BEGIN{OFS=\"\t\"}
       {print $1, $2, $3, \".\", \".\", $4, $5}' unassigned_reads_clusters.bed > unassigned_clusters.bed")
system(paste("bedtools closest -a unassigned_clusters.bed -b gene_ranges_sorted.bed -g ",
             genome_index, " -s -D a -id -fu > clusters_results.txt"))
system("sort -k1,1 -k2,2n intergenic_clusters.bed > intergenic_reads_clusters_sorted.bed")
system(paste("bedtools closest -a intergenic_reads_clusters_sorted.bed -b gene_ranges_sorted.bed -s -D a -id -fu > clusters_results.txt"))

clusters_data <- read.table("clusters_results.txt", sep = "\t")

clusters_data <- clusters_data[clusters_data$V17 != '.', ]

# Sort by cluster sizes:
clusters_data <- clusters_data[order(clusters_data$V7, decreasing = TRUE), ]

# Leave only informative columns:
clusters_data <- clusters_data[, c(1,2,3,6,7,17,18)]

# Rename columns:
colnames(clusters_data) <- c("chr","start","end","strand","count","closest_gene","distance")

rownames(clusters_data) <- seq(1, nrow(clusters_data))

print("Largest intergenic clusters:")
print(head(clusters_data))

# Plot histogram of clusters (of size bigger than 10)
distances <- clusters_data$distance[clusters_data$count > 10] * -1
distances <- sort(distances)
hist(distances[distances < 100000], main = "Intergenic Read Clusters",
     xlab = "Distance from 3' ends", ylab = "Count of clusters")

# Save lists:
write.csv(clusters_data, "intergenic_clusters.csv")

# Clean up:
system("rm unassigned_reads_clusters.bed unassigned_clusters_filtered.bed")
system("rm clusters_results.txt")
system("rm intergenic_clusters.bed")
system("rm intergenic_reads_clusters.bed")
system("rm intergenic_reads_clusters_sorted.bed")

print("Done.")
