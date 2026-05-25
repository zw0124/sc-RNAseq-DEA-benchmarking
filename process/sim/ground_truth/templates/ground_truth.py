#!/usr/bin/env python3

import scanpy as sc


file_in = "${input_anndata}"
file_out = "de_${meta.scenario}_${meta.run}.tsv"
scenario = "${meta.scenario}"

print(f'Input: {file_in}')
print(f'Output: {file_out}')
print(f'Scenario: {scenario}')

adata = sc.read_h5ad(file_in)

adata.var[['ConditionDE.Condition1', 'ConditionDE.Condition2']].to_csv(file_out, sep = '\\t')

