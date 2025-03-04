import re
import sys

def find_tata_motif(sequence, strand):
    if strand == "+":
        seq = sequence[:1000]
        pattern = r'TATA[AT]A[AT]'  # Regex for TATA(A/T)A(A/T)
        matches = [(m.start(), m.group()) for m in re.finditer(pattern, seq)]
        closest_match = max(matches, key=lambda x: x[0]) if matches else None
        return 1000 - closest_match[0]
    else:
        seq = sequence[-1000:]
        pattern = r'[AT]T[AT]TATA'  # Regex for TATA(A/T)A(A/T)
        matches = [(m.start(), m.group()) for m in re.finditer(pattern, seq)]
        closest_match = min(matches, key=lambda x: x[0]) if matches else None
        return closest_match

def find_ATrich_segment_downstream(sequence, strand, min_length=10):
    if strand == "+":
        seq = sequence[-1050:]
        for i in range(len(seq) - min_length + 1):
            window = seq[i:i + min_length]
            if window.count('A') >= min_length - 1:  # Allow one mismatch
                match = re.match(r"A+[^A]?A*", seq[i:])  # Extend beyond min_length
                if match:
                    return i - 50, len(match.group())
        return ".", "."

    else:
        seq = sequence[:1050]
        for i in range(len(seq) - min_length, -1, -1):
            window = seq[i:i + min_length]
            if window.count('T') >= min_length - 1:  # Allow one mismatch
                match = re.match(r"T*[^T]?T+", seq[i:])  # Extend beyond min_length
                if match:
                    return 1000 - i - len(match.group()), len(match.group())
        return ".", "."

# file contains sequences of regions + 1000bp around regions
sequence_file = sys.argv[1]
output_file = sys.argv[2]

with open(sequence_file, "r") as sequences, open(output_file, "w") as output:
    for line in sequences:
        region, strand, sequence = line.strip().split("\t")

        tata_distance = find_tata_motif(sequence, strand)
        AT_distance_downstream, segment_len = find_ATrich_segment_downstream(sequence, strand)

        output.write(f"{region}\t{tata_distance}\t{AT_distance_downstream}\t{segment_len}\n")


