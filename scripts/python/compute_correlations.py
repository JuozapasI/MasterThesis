import scanpy as sc
import scipy.stats as stats
import numpy as np

intergenic_list_path = '../../data/downstream/intergenic/filtered_for_correlations.bed'

samples = [
    "PBMC_10x", "PBMC_10x_2", "PBMC_10x_3", "PBMC_indrops", "PBMC_indrops_2", 
    "brain", "brain_2", "eye", "eye_2", "eye_3", "lung_2", "lung_5", "lung_7", "lung_8"
]

adatas = {}
print("Loading samples")
for sample in samples:
    adatas[sample] = sc.read_h5ad(f'../../data/downstream/adatas/{sample}.h5ad')
    adatas[sample] = adatas[sample].raw.to_adata()
    sc.pp.normalize_total(adatas[sample], target_sum=None)
    sc.pp.log1p(adatas[sample])

adatas_intergenic = {}
print("Loading intergenic samples")
for sample in samples:
    adatas_intergenic[sample] = sc.read_h5ad(f'../../data/downstream/adatas/intergenic_combined/{sample}.h5ad')

print("Adjusting adatas")
for sample in samples:
    # Ensure that both adatas contain same cells and same order of them
    mask =  [True if i in adatas[sample].obs_names.values else False for i in adatas_intergenic[sample].obs_names.values]
    adatas_intergenic[sample] = adatas_intergenic[sample][mask]
    adatas[sample] = adatas[sample][adatas_intergenic[sample].obs_names]

print("Computing correlations")
with open(intergenic_list_path, "r") as intergenic_list, open(f'{intergenic_list_path.removesuffix(".bed")}_with_correlations.bed', "w") as output:
    for line in intergenic_list:
        entries = line.strip().split('\t')
        # Intergenic gene name
        gene0 = entries[3]
        # Gene names of closest on the same strand and closest on the opposite strand
        gene1 = entries[9].split(":")[1]
        gene2 = entries[11].split(":")[1]

        r1, r2= {}, {}

        for sample in samples:
            if gene1 in adatas[sample].var_names:
                expr1 = adatas_intergenic[sample][:, gene0].X.toarray().flatten()
                expr2 = adatas[sample][:, gene1].X.toarray().flatten()
                # Check if expression is non constant (particularly 0), if it is, then correlation is not defined
                if not (np.all(expr1 == expr1[0]) or np.all(expr2 == expr2[0])):
                    r, p = stats.spearmanr(expr1, expr2)
                    # Focus only on reliable results
                    if p < 0.05:
                        r1[sample] = r

            if gene2 in adatas[sample].var_names:
                expr1 = adatas_intergenic[sample][:, gene0].X.toarray().flatten()
                expr2 = adatas[sample][:, gene2].X.toarray().flatten()
                # Check if expression is non constant (particularly 0), if it is, then correlation is not defined
                if not (np.all(expr1 == expr1[0]) or np.all(expr2 == expr2[0])):
                    r, p = stats.spearmanr(expr1, expr2)
                    # Focus only on reliable results
                    if p < 0.05:
                        r2[sample] = r

        if(len(r1) == 0):
            max_sample_1 = "."
            max_value_1 = "."
        else:
            max_sample_1 = max(r1, key=r1.get)
            max_value_1 = r1[max_sample_1]

        if(len(r2) == 0):
            max_sample_2 = "."
            max_value_2 = "."
        else:
            max_sample_2 = max(r2, key=r2.get)
            max_value_2 = r2[max_sample_2]

        output.write("\t".join(entries) + f"\t{max_sample_1}\t{max_value_1}\t{max_sample_2}\t{max_value_2}\n")

print("Done!")

