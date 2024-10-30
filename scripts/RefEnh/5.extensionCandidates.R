print("5. Generating list of gene extension candidates.")

# Load intergenic clusters data:
intergenic_clusters <- read.csv("intergenic_clusters.csv", row.names = 1)

# Find extension candidates based on found intergenic clusters:
extension_candidates = intergenic_clusters[intergenic_clusters$gene != ".",]

extension_candidates$distance = -1*extension_candidates$distance

# Leave only gene with clusters sufficiently close:
extension_candidates = extension_candidates[extension_candidates$distance < extensionThreshold, ]

# Compute counts per millions score
extension_candidates$cpm = extension_candidates$count/seq_depth*1000000

# Filter based on cpm:
extension_candidates = extension_candidates[extension_candidates$cpm > cpmThreshold, ]

# Make sure that we have unique entries, in the case there are 2 clusters for one gene, leave only the bigger one).
extension_candidates = extension_candidates[!duplicated(extension_candidates$gene), ]

# Save the list:
write.csv(extension_candidates, "extension_candidates.csv")

# Additionally let's create a list of new gene candidates:
new_genes <- intergenic_clusters[-1*intergenic_clusters$distance > extensionThreshold, ]

# Save the list:
write.csv(new_genes, "new_gene_candidates.csv")

print("Done.")
