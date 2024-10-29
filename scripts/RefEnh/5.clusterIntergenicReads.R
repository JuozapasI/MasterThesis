print("5. Clustering intergenic reads.")

# Get clusters of intergenic reads:
system("bedtools merge -s -c 6 -o distinct,count < intergenic_reads_sorted.bed > intergenic_reads_clusters.bed")
# Add empty columns to make it compatible with bedtools:
system("awk 'BEGIN{OFS=\"\t\"}
       {print $1, $2, $3, \".\", \".\", $4, $5}' intergenic_reads_clusters.bed > intergenic_clusters.bed")
system("sort -k1,1 -k2,2n intergenic_clusters.bed > intergenic_reads_clusters_sorted.bed")
system(paste("bedtools closest -a intergenic_reads_clusters_sorted.bed -b gene_ranges_sorted.bed -g ",
             genome_index, " -s -D a -id -fu > clusters_results.txt"))

clusters_data <- read.table("clusters_results.txt", sep = "\t")

unassigned_clusters <- clusters_data[clusters_data$V17 == '.', ]
clusters_data <- clusters_data[clusters_data$V17 != '.', ]

# Sort by cluster sizes:
unassigned_clusters <- unassigned_clusters[order(unassigned_clusters$V7, decreasing = TRUE), ]
clusters_data <- clusters_data[order(clusters_data$V7, decreasing = TRUE), ]

rownames(clusters_data) <- seq(1, nrow(clusters_data))
rownames(unassigned_clusters) <- seq(1, nrow(unassigned_clusters))

print("Largest intergenic clusters:")
print(head(clusters_data[,c("V1", "V2", "V3", "V6", "V7", "V17")]))

# Plot histogram of clusters (of size bigger than 10)
distances <- clusters_data$V18[clusters_data$V7 > 10] * -1
distances <- sort(distances)
hist(distances[distances < 100000], main = "Intergenic Read Clusters",
     xlab = "Distance from 3' ends", ylab = "Count of clusters")

# Save lists:
write.csv(clusters_data, "intergenic_clusters.csv")
write.csv(unassigned_clusters, "unassigned_clusters.csv")

# Clean up:
system("rm clusters_results.txt")
system("rm intergenic_clusters.bed")
system("rm intergenic_reads_clusters.bed")
system("rm intergenic_reads_clusters_sorted.bed")

print("Done.")