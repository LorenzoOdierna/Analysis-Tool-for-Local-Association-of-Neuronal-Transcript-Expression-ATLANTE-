import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from numba import jit
import warnings

CONFIG = {
    'working_directory': 'C:/Users/godierna/Datafiles',
    'rna_brain_file': 'rna_brain_HPA.tsv',
    'regional_expression_file': 'HPA_RegionalGeneExpression_2025.csv',
    'target_genes_file': 'target_genes.csv',
    'target_genes_column': 'GENES',
    'percentile_threshold': 95,
    'n_target_genes':,  # input lenth of gene list here
'n_simulations': 1000000,
'chunk_size': 1000,
'confidence_multiplier': 3.471,
'output_filtered_genes': 'filtered_target_genes.csv',
'output_frequencies': 'observed_frequencies.csv',
'output_final_results': 'enrichment_analysis_results.csv',
}

os.chdir(CONFIG['working_directory'])

Data_raw = pd.read_table(CONFIG['rna_brain_file'])
nTPM = Data_raw.drop(['Gene', 'TPM', 'pTPM'], axis=1)
BrainRegions_temp = nTPM['Subregion'].unique()
all_genes = nTPM['Gene name'].unique()
filtered_data = []
for gene in tqdm(all_genes, desc="Filtering genes"):
    gene_data = nTPM[nTPM['Gene name'] == gene].copy()
    threshold = np.percentile(gene_data['nTPM'], CONFIG['percentile_threshold'])
    gene_filtered = gene_data[gene_data['nTPM'] >= threshold]
    filtered_data.append(gene_filtered)
Result_preprocessed = pd.concat(filtered_data, ignore_index=True)
Result_preprocessed.to_csv(CONFIG['regional_expression_file'], index=False)


@jit(nopython=True)
def count_frequencies_numba(subregion_indices, brain_region_indices, n_regions):
    counts = np.zeros(n_regions)
    for idx in range(len(subregion_indices)):
        region_idx = brain_region_indices[idx]
        if 0 <= region_idx < n_regions:
            counts[region_idx] += 1
    return counts


Data = pd.read_table(CONFIG['rna_brain_file'])
Result = pd.read_csv(CONFIG['regional_expression_file'])
target_genes = pd.read_csv(CONFIG['target_genes_file'])

BrainRegions = Result['Subregion'].unique()
region_to_idx = {region: idx for idx, region in enumerate(BrainRegions)}

Target_Genes = Result[Result['Gene name'].isin(target_genes[CONFIG['target_genes_column']])]
Target_Genes.to_csv(CONFIG['output_filtered_genes'])

n_genes = CONFIG['n_target_genes']

target_subregion_indices = np.array([region_to_idx.get(region, -1) for region in Target_Genes['Subregion']])
brain_region_indices = np.array([region_to_idx[region] for region in Target_Genes['Subregion']])
numbers = count_frequencies_numba(target_subregion_indices, brain_region_indices, len(BrainRegions))

freqdata = pd.DataFrame({
    'Frequency': numbers,
    'BrainRegion': BrainRegions
})
freqdata.to_csv(CONFIG['output_frequencies'])

humangenome = Data['Gene name'].unique()
result_genes = Result['Gene name'].unique()
gene_to_idx = {gene: idx for idx, gene in enumerate(result_genes)}

gene_subregions = {}
for gene in tqdm(humangenome, desc="Building gene mappings"):
    gene_rows = Result[Result['Gene name'] == gene]
    if not gene_rows.empty:
        subregions = gene_rows['Subregion'].values
        gene_subregions[gene] = [region_to_idx.get(sr, -1) for sr in subregions]

n_sets = CONFIG['n_simulations']
chunk_size = CONFIG['chunk_size']
results = np.zeros((len(BrainRegions), n_sets))

for chunk_start in tqdm(range(0, n_sets, chunk_size), desc="Processing chunks"):
    chunk_end = min(chunk_start + chunk_size, n_sets)
    chunk_size_actual = chunk_end - chunk_start

    for i in range(chunk_size_actual):
        random_genes = np.random.choice(humangenome, size=n_genes, replace=False)
        chunk_counts = np.zeros(len(BrainRegions))

        for gene in random_genes:
            if gene in gene_subregions:
                for region_idx in gene_subregions[gene]:
                    if region_idx >= 0:
                        chunk_counts[region_idx] += 1

        results[:, chunk_start + i] = chunk_counts

means = np.mean(results, axis=1)
stdevs = np.std(results, axis=1)

output = pd.DataFrame({
    'BrainRegion': freqdata['BrainRegion'],
    'Means': means,
    'Stdevs': stdevs,
    'MOE': stdevs * CONFIG['confidence_multiplier'],
})

output['Upper'] = dff['Means'] + dff['MOE']
output['Lower'] = dff['Means'] - dff['MOE']
output['Observed'] = freqdata['Frequency']
output['SigHigh'] = dff['Observed'] > dff['Upper']

dff.to_csv(CONFIG['output_final_results'])