import sys

# Check for correct number of arguments
if len(sys.argv) != 3:
    print("Usage: python convert_ncbi_to_gencode.py <input_gtf> <output_gtf>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Mapping NCBI RefSeq chromosome names to GENCODE chromosome names
ncbi_to_gencode = {
    "NC_000001.11": "1",
    "NC_000002.12": "2",
    "NC_000003.12": "3",
    "NC_000004.12": "4",
    "NC_000005.10": "5",
    "NC_000006.12": "6",
    "NC_000007.14": "7",
    "NC_000008.11": "8",
    "NC_000009.12": "9",
    "NC_000010.11": "10",
    "NC_000011.10": "11",
    "NC_000012.12": "12",
    "NC_000013.11": "13",
    "NC_000014.9": "14",
    "NC_000015.10": "15",
    "NC_000016.10": "16",
    "NC_000017.11": "17",
    "NC_000018.10": "18",
    "NC_000019.10": "19",
    "NC_000020.11": "20",
    "NC_000021.9": "21",
    "NC_000022.11": "22",
    "NC_000023.11": "X",
    "NC_000024.10": "Y",
    "NC_012920.1": "M"  # Mitochondrial chromosome
}

# Open the input and output files
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            # Write header lines as-is
            outfile.write(line)
        else:
            # Process annotation lines
            fields = line.strip().split("\t")
            chrom = fields[0]
            # Replace chromosome name if in the mapping
            if chrom in ncbi_to_gencode:
                fields[0] = ncbi_to_gencode[chrom]
                # Write the modified line to the output
                outfile.write("\t".join(fields) + "\n")

print(f"Chromosome names converted and saved to {output_file}")

