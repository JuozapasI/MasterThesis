import pandas as pd
import sys

def load_bedgraph(file_path):
    """Load a BEDGraph file into a pandas DataFrame."""
    bedgraph = pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end", "coverage"])
    bedgraph["chrom"] = bedgraph["chrom"].astype(str)
    return bedgraph

def load_bed(file_path):
    """Load a BED file into a pandas DataFrame."""
    bed = pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end", "strand", "count"])
    bed["chrom"] = bed["chrom"].astype(str)
    return bed

def find_subregion(bedgraph, bed_region, strand):
    """Find subregion with relatively high coverage for a single region."""
    chrom, start, end = bed_region["chrom"], bed_region["start"], bed_region["end"]
    
    # Filter BEDGraph data for the relevant region
    subgraph = bedgraph[(bedgraph["chrom"] == chrom) & (bedgraph["start"] < end) & (bedgraph["end"] > start)]
    
    if subgraph.empty:
        return -1, -1  # No data in the region

    # Find the bin with maximum coverage
    max_row = subgraph.loc[subgraph["coverage"].idxmax()]
    max_coverage = max_row["coverage"]
    max_start, max_end = max_row["start"], max_row["end"]

    # Expand to adjacent regions with coverage >= max_coverage / 2
    threshold = max_coverage / 2

    # Expand left
    left_start = max_start
    for _, row in subgraph[subgraph["end"] <= max_start].iloc[::-1].iterrows():
        if row["coverage"] >= threshold:
            left_start = row["start"]
        else:
            break

    # Expand right
    right_end = max_end
    for _, row in subgraph[subgraph["start"] >= max_end].iterrows():
        if row["coverage"] >= threshold:
            right_end = row["end"]
        else:
            break

    return left_start, right_end

def process_regions(bed_file, forward_bedgraph, reverse_bedgraph, output_file):
    """Process all regions in the BED file and output refined subregions."""
    regions = load_bed(bed_file)
    forward = load_bedgraph(forward_bedgraph)
    reverse = load_bedgraph(reverse_bedgraph)

    results = []

    for _, region in regions.iterrows():
        strand = region["strand"]
        if strand == "+":
            left, right = find_subregion(forward, region, strand)
        elif strand == "-":
            left, right = find_subregion(reverse, region, strand)
        else:
            continue
        results.append([region["chrom"], left, right, region["strand"], region["count"]])
    
    results = pd.DataFrame(results, columns=["chrom", "start", "end", "strand", "id"])
    results["name"] = "."
    results["score"] = "."
    # Sort as bed file ('sort -k1,1 -k2,2n')
    sorted_df = df.sort_values(
    by=['chrom', 'start'],
    key=lambda col: col.astype(str) if col.name == 'chrom' else col,
    ascending=[True, True]
    )
    # Write refined regions to output file
    results.to_csv(output_file, sep="\t", header=False, index=False, columns=["chrom", "start", "end", "name", "score", "strand", "id"])

if __name__ == "__main__":
    # Get command-line arguments
    bed_file = sys.argv[1]
    forward_bedgraph = sys.argv[2]
    reverse_bedgraph = sys.argv[3]
    output_file = sys.argv[4]

    # Process regions
    process_regions(bed_file, forward_bedgraph, reverse_bedgraph, output_file)

