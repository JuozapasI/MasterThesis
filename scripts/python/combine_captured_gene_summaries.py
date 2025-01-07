import pandas as pd
import os
import re
import sys

# Directory where the CSV files are stored
dir_path = sys.argv[1]
output = sys.argv[2]

# Automatically detect sample names based on the files in the directory
sample_files = [f for f in os.listdir(dir_path) if f.endswith(".10x.csv")]
samples = sorted(set(re.sub(r'\.10x\.csv$', '', f) for f in sample_files))

# Initialize an empty DataFrame to store the combined data
gene_types = set()
combined_data = {}

# Read each sample's files
for sample in samples:
    # Initialize a dictionary to store gene data for this sample
    sample_data = {}

    # Read the .10x.csv file
    with open(os.path.join(dir_path, f"{sample}.10x.csv"), "r") as f:
        for line in f:
            parts = line.strip().split()
            score, gene_type = parts[0], parts[1]
            gene_types.add(gene_type)
            sample_data[f"{gene_type}.10x"] = score

    # Read the .final.csv file
    with open(os.path.join(dir_path, f"{sample}.final.csv"), "r") as f:
        for line in f:
            parts = line.strip().split()
            score, gene_type = parts[0], parts[1]
            gene_types.add(gene_type)
            sample_data[f"{gene_type}.final"] = score

    # Store the data for this sample
    combined_data[sample] = sample_data

# Prepare the LaTeX table header
header = ["gene_type"] + [f"{sample}" for sample in samples]

# Prepare the data rows for the LaTeX table
table_rows = []

# Sort the gene types alphabetically
for gene_type in sorted(gene_types):
    row = [gene_type]
    
    # For each sample, check if the gene_type exists and add the appropriate score
    for sample in samples:
        score_10x = combined_data[sample].get(f"{gene_type}.10x", "0")
        score_final = combined_data[sample].get(f"{gene_type}.final", "0")
        
        # Combine the scores in the required format
        row.append(f"{score_10x} ({score_final})")
    
    table_rows.append(row)

latex_table = pd.DataFrame(table_rows, columns=header)
latex_table['SortKey'] = latex_table.iloc[:, 1].str.extract(r'(\d+)').astype(int)
latex_table = latex_table.sort_values(by='SortKey', ascending=False)
latex_table = latex_table.drop(columns=['SortKey'])


# Output the LaTeX table
latex_code = latex_table.to_latex(index=False, escape=False)
latex_code = latex_code.replace("\\begin{tabular}", "\\begin{table}[htbp]\n\\centering\n\\begin{tabular}")
latex_code = latex_code.replace("\\end{tabular}", "\\end{tabular}\n\\caption{Combined captured gene types summary}\n\\label{tab:geneTypes}\n\\end{table}")
latex_code = latex_code.replace("_", r"\_")


with open(output, "w") as output_file:
    output_file.write(latex_code)

