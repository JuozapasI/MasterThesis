import os
import pandas as pd
import sys

def read_file(filepath):
    """Read a file and return a dictionary of category counts."""
    data = {}
    with open(filepath, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) >= 3:
                count, source, category = int(parts[0]), parts[1], " ".join(parts[2:])
                key = f"{source} {category}"
                data[key] = count
    return data

def generate_latex_table(directory):
    """Generate a LaTeX table from files in the given directory."""
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and f.endswith(".csv")]
    
    # Remove suffixes from filenames
    sample_names = [os.path.splitext(f)[0] for f in files]

    # Read all files and combine data
    combined_data = {}
    for filename in files:
        filepath = os.path.join(directory, filename)
        file_data = read_file(filepath)
        combined_data[filename] = file_data

    # Get all unique categories
    all_categories = sorted(set(key for file_data in combined_data.values() for key in file_data.keys()))
    
    # Create a DataFrame
    df = pd.DataFrame(index=all_categories, columns=sample_names).fillna(0)
    for filename, file_data in combined_data.items():
        sample_name = os.path.splitext(filename)[0]
        for category, count in file_data.items():
            df.at[category, sample_name] = count

    # Split categories into primary and subcategories
    df['Primary'] = [cat.split()[0] for cat in df.index]
    df['Category'] = [cat[len(primary) + 1:] for cat, primary in zip(df.index, df['Primary'])]
    
    # Organize data for LaTeX
    df = df.sort_values(by=['Primary', df.columns[0]], ascending=[True,False])

    latex_data = []
    for primary, group in df.groupby('Primary', sort=False):
        latex_data.append([primary, *["" for _ in sample_names]])  # Primary row
        for _, row in group.iterrows():
            latex_data.append([row['Category'], *[int(row[file]) for file in sample_names]])
    
    # Generate LaTeX table
    latex_table = "\\begin{table}[htbp]\n\\centering\n\\small\n\\begin{tabular}{l" + "r" * len(sample_names) + "}\n"
    latex_table += "Sample & " + " & ".join([f"\\texttt{{{name}}}" for name in sample_names]) + " \\\\\n\\hline\n"
    for row in latex_data:
        latex_table += " & ".join(map(str, row)) + " \\\\\n"
    latex_table += "\\end{tabular}\n\\caption{Summary of gene types intersecting with unassigned reads}\n\\labe{tab:intersectingGenes}\n\\end{table}"

    latex_table = latex_table.replace("_", r"\_")
    
    return latex_table


directory = sys.argv[1]
output = sys.argv[2]

# Generate LaTeX table
latex_table = generate_latex_table(directory)

# Save the LaTeX table to a file
with open(output, "w") as output_file:
    output_file.write(latex_table)

