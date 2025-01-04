import os
import pandas as pd
import sys

def process_summary_file(file_path):
    data = []
    data.append(["sample", os.path.splitext(os.path.basename(file_path))[0]])
    with open(file_path, 'r') as f:
        lines = f.readlines()
        # Extract total reads
        for line in lines:
            if not line.strip():
                continue
            parts = line.strip().split(',')
            if "Total reads" in line:
                data.append(["Total Reads", int(parts[1])])
            elif ":" in parts[0]:
                label = parts[0].replace(':', '').strip()
                percentage = round(float(parts[2]), 2) if len(parts) > 2 else 0
                data.append([label, percentage])
            else:
                data.append([parts[0], ""])
    return data

def generate_latex_table(data, output_file):
    latex_table = "\\begin{table}[htbp]\n\\centering\n\\begin{tabular}{|l|" + "c|" * (len(data[0]) - 1) + "}\n"
    latex_table += "\\hline\n"
    for row in data:
        latex_table += " & ".join([str(i) for i in row]) + " \\\\\n"
    
    latex_table += "\\hline\n\\end{tabular}\n\\caption{Summaries of Reads and Percentages.}\n\\label{tab:countSummariesCombined}\n\\end{table}"
    
    with open(output_file, "w") as f:
        f.write(latex_table)

# Main logic
def main(input_dir, output_file):
    combined_data = []
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".csv"):  # Assuming summary files are .csv
            file_path = os.path.join(input_dir, file_name)
            file_data = process_summary_file(file_path)
            combined_data.append(file_data)
    
    data = combined_data[0]
    for i in range(1, len(combined_data)):
        for j in range(len(combined_data[0])):
            data[j].append(combined_data[i][j][1])

    
    generate_latex_table(data, output_file)

input_directory = sys.argv[1]
output_latex_file = sys.argv[2]
main(input_directory, output_latex_file)



