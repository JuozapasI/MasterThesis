print("3. Isolating intergenic reads.")

# Firstly, let's take those reads that were not assigned to any genes (i.e. has "GN:Z:-" tag).
# Additionally, let's filter those that have corrupted barcodes or UMIs.
# Finally, filter duplicates based on barcode and UMI combination
# -F 256 and -F 1024 filter duplicate and secondary allignment reads
system(paste("( samtools view -H ", bam_file, " && samtools view -F 256 -F 1024 ", bam_file,
    " | grep \"GN:Z:-\" | grep -E \"UB:Z:[A-Z]{", umi_length, "}\" | grep -E \"CB:Z:[A-Z]{", barc_length,
    "}\" | awk 'BEGIN { FS=\"\t\"; OFS=\"\t\" } { ", 
    "if ($0 ~ /CB:Z:/ && $0 ~ /UB:Z:/) { ",
    "split($0, a, \"CB:Z:\"); split(a[2], b, \" \"); cb = b[1]; ",
    "split($0, a, \"UB:Z:\"); split(a[2], b, \" \"); ub = b[1]; ",
    "if (!seen[cb\":\"ub]++) { print; } } }' ; ) | samtools view -h -b - > unassigned_reads.bam",
    sep = ""))

# Convert to bed
system("bedtools bamtobed -i unassigned_reads.bam > unassigned_reads.bed")
# Sort
system("sort -k1,1 -k2,2n unassigned_reads.bed > unassigned_reads_sorted.bed")
# Take only intergenic (not intersecting with genes) reads
system("bedtools intersect -v -a unassigned_reads.sorted.bed -b gene_ranges_sorted.bed > intergenic_reads.bed")
# Sort
system("sort -k1,1 -k2,2n intergenic_reads.bed > intergenic_reads_sorted.bed")
# Save also the unassigned reads that interesct with genes:
system("bedtools intersect -a unassigned_reads.sorted.bed -b gene_ranges_sorted.bed > intersecting_reads.bed")
# Sort
system("sort -k1,1 -k2,2n intersecting_reads.bed > intersecting_reads_sorted.bed")

# Clean-up
system("rm intersecting_reads.bed intergenic_reads.bed unassigned_reads.bed unassigned_reads_sorted.bed")

print("Done.")
