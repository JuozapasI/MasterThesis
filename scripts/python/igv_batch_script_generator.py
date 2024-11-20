import sys
import os

# Check if the correct number of arguments is provided
if len(sys.argv) != 8:
    print("Usage: python generate_igv_batch_script.py <bed_file> <bam_file> <output_folder> <genome_file> <genome_name> <annotation_file> <batch_script>")
    sys.exit(1)

# Get the arguments from sys.argv
bed_file = sys.argv[1]  # Input BED file
bam_file = sys.argv[2]  # BAM file path
output_folder = sys.argv[3]  # Directory to save the snapshot images
genome_file = sys.argv[4]  # Path to custom genome FASTA file
genome_name = sys.argv[5]  # Name of the genome (e.g., "custom_genome")
annotation_file = sys.argv[6]  # Path to the transcriptomic annotation file (GTF or GFF)
batch_script = sys.argv[7]  # Path to output the IGV batch script

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Generate the batch script
with open(bed_file, "r") as bed, open(batch_script, "w") as script:
    script.write(f"new\n")
    script.write(f"genome {genome_file} {genome_name}\n")  # Load the custom genome
    script.write(f"load {bam_file}\n")  # Load the BAM file
    script.write(f"load {annotation_file}\n")  # Load the transcriptomic annotation file
    script.write(f"snapshotDirectory {output_folder}\n\n")
    
    # Write commands for each region in the BED file
    for line in bed:
        chrom, start, end = line.strip().split()[:3]
        script.write(f"goto {chrom}:{int(start) - 1000}-{int(end) + 1000}\n")
        script.write(f"snapshot {chrom}_{start}_{end}.png\n\n")

print(f"Batch script saved to {batch_script}")
