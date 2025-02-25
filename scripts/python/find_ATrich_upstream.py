import sys

def compute_score(seq):
    """Computes the sliding window score for 10-mers, allowing 1 mistake."""
    for i in range(len(seq) - 9):
        window = seq[i:i+10]
        if window.count('A') >= 9 or window.count('T') >= 9:
            return i
    return -1


def process_files(tsv_file, bed_file, output_file):
    # Read sequences from TSV
    seq_dict = {}
    with open(tsv_file, 'r') as tsv:
        for line in tsv:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                seq_dict[parts[0]] = parts[1]

    # Process BED file and compute scores
    with open(bed_file, 'r') as bed, open(output_file, 'w') as out:
        for line in bed:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            chrom, start, end, name, _, strand = parts[:6]
            start, end = int(start), int(end)
            
            if name not in seq_dict:
                continue
            
            seq = seq_dict[name]
            score_pos = compute_score(seq[::-1] if strand == '+' else seq)
            
            out.write(f"{name}\t{score_pos}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <tsv_file> <bed_file> <output_file>")
        sys.exit(1)
    
    tsv_file = sys.argv[1]
    bed_file = sys.argv[2]
    output_file = sys.argv[3]
    
    process_files(tsv_file, bed_file, output_file)