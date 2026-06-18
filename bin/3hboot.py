#!/usr/bin/env python3

import argparse

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from tqdm import tqdm

np.random.seed(0)


parser = argparse.ArgumentParser(
    prog="Hierarchical Bootstrap",
    description="Original recursive hierarchical bootstrap with cell-level resampling"
)
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

file_in = args.input
file_out = args.output

GROUP_COL = "Condition"
HIERARCHY_PRIORITY = ["Batch", "Sample"]
BOOTSTRAP_ITER = 10_000


def detect_hierarchy(adata):
    hierarchy = []
    for col in HIERARCHY_PRIORITY:
        if col not in adata.obs.columns:
            continue
        if col == "Batch" and adata.obs[col].astype(str).nunique() <= 1:
            print("[INFO] Skipping Batch hierarchy level (missing variation: unique<=1).")
            continue
        hierarchy.append(col)

    print("[INFO] Detected hierarchy:", hierarchy)
    return hierarchy


def sample_indices(adata, mask, hierarchy):
    if len(hierarchy) == 0:
        cells = np.where(mask)[0]
        return np.random.choice(cells, size=len(cells), replace=True)

    level = hierarchy[0]
    units = adata.obs.loc[mask, level].unique()
    sampled_units = np.random.choice(units, size=len(units), replace=True)

    idx = []
    for unit in sampled_units:
        next_mask = ((adata.obs[level] == unit).to_numpy() & mask)
        idx.extend(sample_indices(adata, next_mask, hierarchy[1:]))

    return idx


def row_mean(X, idx):
    mean_expr = X[idx, :].mean(axis=0)
    return np.asarray(mean_expr).ravel()


def hierarchical_bootstrap(adata, groupby, hierarchy, n_iter):
    X = adata.X
    pseudo_samples = []
    pseudo_obs = []
    groups = adata.obs[groupby].unique()

    for i in tqdm(range(n_iter), desc="Hierarchical bootstrap"):
        for group in groups:
            mask = (adata.obs[groupby] == group).to_numpy()
            idx = sample_indices(adata, mask, hierarchy)

            pseudo_samples.append(row_mean(X, idx))
            pseudo_obs.append({
                "bootstrap_id": f"{group}_{i}",
                groupby: group
            })

    return AnnData(
        X=np.vstack(pseudo_samples),
        var=adata.var.copy(),
        obs=pd.DataFrame(pseudo_obs).set_index("bootstrap_id")
    )


def bootstrap_pvalue(adata, gene, groupby):
    groups = adata.obs[groupby].unique()
    x = np.asarray(adata[:, gene].X).ravel()

    x1 = x[(adata.obs[groupby] == groups[0]).to_numpy()]
    x2 = x[(adata.obs[groupby] == groups[1]).to_numpy()]

    n = min(len(x1), len(x2))
    wins = np.sum(x1[:n] > x2[:n])

    return min(1.0, 2 * min(wins / n, 1 - wins / n))


def calculate_observed_lfc(adata, groupby):
    groups = adata.obs[groupby].unique()
    cond_a, cond_b = groups[0], groups[1]

    mask_a = (adata.obs[groupby] == cond_a).to_numpy()
    mask_b = (adata.obs[groupby] == cond_b).to_numpy()

    mean_a = row_mean(adata.X, mask_a)
    mean_b = row_mean(adata.X, mask_b)

    return pd.Series(
        np.log2((mean_a + 1e-9) / (mean_b + 1e-9)),
        index=adata.var_names
    )


if __name__ == "__main__":
    print("[INFO] Loading AnnData...")
    adata = sc.read_h5ad(file_in)

    print("[INFO] Applying Total Count Normalization...")
    sc.pp.normalize_total(adata, target_sum=1e6)

    hierarchy = detect_hierarchy(adata)
    observed_lfc = calculate_observed_lfc(adata, GROUP_COL)

    print(f"\n[INFO] Starting original hierarchical bootstrap ({BOOTSTRAP_ITER:,} Iters)...")
    adata_hb = hierarchical_bootstrap(
        adata,
        GROUP_COL,
        hierarchy,
        BOOTSTRAP_ITER
    )

    pvals = {}
    for gene in tqdm(adata.var_names, desc="Computing p-values"):
        pvals[gene] = bootstrap_pvalue(adata_hb, gene, GROUP_COL)

    pvals = pd.Series(pvals)

    pd.DataFrame({
        "gene": pvals.index,
        "p_val": pvals.values,
        "lfc": observed_lfc.loc[pvals.index].values
    }).to_csv(file_out, sep="\t", index=False)

