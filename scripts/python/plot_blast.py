import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob

def get_tissue_type(sample_name):
    if 'PBMC' in sample_name:
        return 'PBMC'
    elif 'brain' in sample_name:
        return 'brain'
    elif 'lung' in sample_name:
        return 'lung'
    elif 'eye' in sample_name:
        return 'eye'
    return 'Unknown'

id_to_name = {
    ">NC_018476.1": "Simbu virus RdRp",
    ">NC_032111.1": "BeAn 58058 virus",
    ">NC_048651.1": "Aeribacillus phage AP45",
    ">NC_001422.1": "Escherichia phage phiX174",
    ">NC_043329.1": "Diolcogaster facetosa bracovirus",
    ">NC_055235.1": "Papiine betaherpesvirus 3",
    ">NC_008168.1": "Choristoneura fumiferana granulovirus",
    "Other": "Other",
    ">NC_076967.1": "Colobine gammaherpesvirus 1",
    ">NC_022518.1": "Human endogenous retrovirus K113",
    ">NC_028045.1": "Tadarida brasiliensis circovirus 1"
}

# Read and combine all
dfs = []
for fname in glob.glob('/home/juzis/Kita/cluster/data/blast/viral/plot/*_proportions.tsv'):
    df = pd.read_csv(fname, sep='\t', header=None, names=['Sample', 'Virus', 'Proportion'])
    df['Tissue'] = df['Sample'].apply(get_tissue_type)
    df['Virus'] = df['Virus'].map(id_to_name)
    dfs.append(df)
data = pd.concat(dfs)

# Pivot data: rows = samples, columns = viruses
pivot = data.pivot_table(index='Sample', columns='Virus', values='Proportion', fill_value=0)

fig, axes = plt.subplots(1, 4, figsize=(10, 6), width_ratios=[5, 2, 3, 4])

# Tissue types and corresponding sample indices (adjust accordingly)
tissues = ['PBMC', 'Brain', 'Eye', 'Lung']
tissue_samples = {
    'PBMC': [0, 1, 2, 3, 4],  # Adjust indices based on your actual sample order
    'Brain': [5, 6],
    'Eye': [7, 8, 9],
    'Lung': [10, 11, 12, 13]
}

# Iterate over each tissue type and create the same bar plot for each subplot
for i, tissue in enumerate(tissues):
    ax = axes[i]
    
    # Select the relevant samples for the tissue
    tissue_data = pivot.iloc[tissue_samples[tissue]]
    
    # Plot the bar chart in the corresponding subplot
    tissue_data.plot(kind='bar', stacked=True, ax=ax, colormap='Paired', width=0.8)
    
    # Set title for each subplot
    ax.set_title(f'{tissue}')
    
    # Sample names will be taken from input data
    ax.set_xticks(range(len(tissue_samples[tissue])))  # Adjust according to the number of samples
    ax.set_xticklabels(pivot.index[tissue_samples[tissue]], rotation=45)

    # Set y-axis label only on the first subplot
    if i == 0:
        ax.set_ylabel('Proportion')

    # Remove y-axis (ticks, labels, and spine) on all but the first subplot
    if i != 0:
        ax.set_yticks([])  # Remove y-axis ticks
        ax.set_yticklabels([])  # Remove y-axis labels
        ax.spines['left'].set_visible(False)  # Remove y-axis spine

    # Remove borders for all subplots
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if i == 0:
        handles, labels = ax.get_legend_handles_labels()

    ax.legend().set_visible(False)
    ax.set_xlabel(None)

# Add a single legend below all subplots
plt.tight_layout(rect=[0, 0.17, 1, 1])
fig.legend(handles, labels, loc='lower center', ncol=3, title='Virus')

# Save and display the plot
plt.savefig('virus_proportions_by_tissue.png', dpi=300)
plt.show()