#!/usr/bin/env python3

import scanpy as sc
import numpy as np

file_in = "${input_h5ad}"
seed = int("${run}")
sample_size = int("${sample_size}")
file_out = "${meta.scenario}_${meta.run}.h5ad"

print(f"Reading {file_in}...")
adata = sc.read_h5ad(file_in)

if 'Condition' in adata.obs.columns:
    adata = adata[adata.obs['Condition'] == 'normal'].copy()


donors = adata.obs['Sample'].unique()
print(f"Total unique donors: {len(donors)}")

num_donors = sample_size * 2
if len(donors) < num_donors:
    raise ValueError(
        f"Need at least {num_donors} unique donors for {sample_size}v{sample_size} sampling, got {len(donors)}."
    )

print(f"Generating {file_out} with seed={seed}, sample_size={sample_size}v{sample_size}...")

rng = np.random.default_rng(seed)
selected_donors = rng.choice(donors, num_donors, replace=False)
group1 = selected_donors[:sample_size]
group2 = selected_donors[sample_size:]

adata_sub = adata[adata.obs['Sample'].isin(selected_donors)].copy()

adata_sub.obs['Condition'] = 'group1'
adata_sub.obs.loc[adata_sub.obs['Sample'].isin(group2), 'Condition'] = 'group2'

adata_sub.var['ConditionDE.Condition1'] = 1.0
adata_sub.var['ConditionDE.Condition2'] = 1.0

adata_sub.write_h5ad(file_out)
print("Finished!")
