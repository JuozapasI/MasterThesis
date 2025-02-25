import sys

def contains_AT_rich_region(sequence):
    """
    Check if the sequence contains 10 consecutive 'A' or 'T' allowing one mismatch.
    """
    for i in range(len(sequence) - 9):  # Sliding window of size 10
        window = sequence[i:i+10]
        at_count = sum(1 for c in window if c in "AT")
        if at_count >= 9:  # At least 9 out of 10 must be A/T
            return True
    return False

def process_tsv(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            name, sequence = line.strip().split("\t")
            if contains_AT_rich_region(sequence):
                print(f"{name}\t{1}")
            else:
                print(f"{name}\t{0}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input.tsv")
        sys.exit(1)

    process_tsv(sys.argv[1])