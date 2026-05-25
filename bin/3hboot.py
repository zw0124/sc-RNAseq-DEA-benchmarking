#!/usr/bin/env python3

import argparse
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse
from tqdm import tqdm
import gc

np.random.seed(0)

parser = argparse.ArgumentParser(
    prog="three-level hierarchical bootstrap",
    description="Run three-level hierarchical bootstrap using the batch-sample-cell hierarchy."
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
    print(f"[INFO] Detected hierarchy: {hierarchy}")
    return hierarchy

def build_tree_fast(df, cols):
    """
    Convert a dataframe into a batch-sample-cell hierarchy with cell indices at the leaves.
    """
    if len(cols) == 1:
        # The terminal level stores cell indices, which are resampled in sample_tree().
        return {k: v.index.values for k, v in df.groupby(cols[0])}
    else:
        col = cols[0]
        return {k: build_tree_fast(v, cols[1:]) for k, v in df.groupby(col)}

def generate_group_weights_from_tree(tree, n_cells, n_iters):
    """
    Generate bootstrap sampling weights from a precomputed hierarchy.
    """
    W = np.zeros((n_cells, n_iters), dtype=np.int32)
    
    def sample_tree(node):
        if isinstance(node, np.ndarray):
            if len(node) > 0:
                # Cell-level bootstrap: resample cells within the selected sample.
                return np.random.choice(node, size=len(node), replace=True)
            return np.array([], dtype=int)
        
        keys = list(node.keys())
        if len(keys) == 0:
            return np.array([], dtype=int)
            
        sampled_keys = np.random.choice(keys, size=len(keys), replace=True)
        results = [sample_tree(node[k]) for k in sampled_keys]
        
        if results:
            return np.concatenate(results)
        return np.array([], dtype=int)

    for i in range(n_iters):
        sampled_cells = sample_tree(tree)
        if len(sampled_cells) > 0:
            W[:, i] = np.bincount(sampled_cells, minlength=n_cells)
            
    return W

def run_hierarchical_bootstrap(adata, group_col, hierarchy, total_iters, chunk_size=10000, gene_subset_idx=None):
    conds = adata.obs[group_col].unique()
    cond_A, cond_B = conds[0], conds[1]
    
    mask_A = (adata.obs[group_col] == cond_A).values
    mask_B = (adata.obs[group_col] == cond_B).values
    
    X = adata.X.T 
    if gene_subset_idx is not None:
        X = X[gene_subset_idx, :]
        
    if not scipy.sparse.isspmatrix_csr(X):
        X = scipy.sparse.csr_matrix(X)
        
    X_A_T = X[:, mask_A]
    X_B_T = X[:, mask_B]
    
    obs_A = adata.obs.iloc[mask_A].reset_index(drop=True)
    obs_B = adata.obs.iloc[mask_B].reset_index(drop=True)
    
    print("[INFO] Building independent hierarchical trees for Group A and B...")
    tree_A = build_tree_fast(obs_A, hierarchy)
    tree_B = build_tree_fast(obs_B, hierarchy)
    
    total_wins_gt = np.zeros(X.shape[0], dtype=np.int64)
    total_wins_lt = np.zeros(X.shape[0], dtype=np.int64)
    
    n_chunks = int(np.ceil(total_iters / chunk_size))
    
    for i in tqdm(range(n_chunks), desc=f"Bootstrap Chunks (Total: {total_iters})"):
        current_chunk_size = min(chunk_size, total_iters - i * chunk_size)
        
        W_A = generate_group_weights_from_tree(tree_A, len(obs_A), current_chunk_size)
        W_B = generate_group_weights_from_tree(tree_B, len(obs_B), current_chunk_size)
        
        Sum_A = X_A_T.dot(W_A)
        Sum_B = X_B_T.dot(W_B)
        
        N_A = W_A.sum(axis=0)
        N_B = W_B.sum(axis=0)
        N_A[N_A == 0] = 1 
        N_B[N_B == 0] = 1
        
        Mean_A = Sum_A / N_A
        Mean_B = Sum_B / N_B
        
        total_wins_gt += np.sum(Mean_A > Mean_B, axis=1)
        total_wins_lt += np.sum(Mean_A < Mean_B, axis=1)
        
        del W_A, W_B, Sum_A, Sum_B, Mean_A, Mean_B
        gc.collect()
        
    pvals = np.minimum(total_wins_gt, total_wins_lt) * 2.0 / total_iters
    return pvals

def calculate_observed_lfc(adata, group_col):
    conds = adata.obs[group_col].unique()
    cond_A, cond_B = conds[0], conds[1]

    mask_A = (adata.obs[group_col] == cond_A).values
    mask_B = (adata.obs[group_col] == cond_B).values

    X = adata.X
    if not scipy.sparse.issparse(X):
        X = scipy.sparse.csr_matrix(X)

    mean_A = np.asarray(X[mask_A, :].mean(axis=0)).ravel()
    mean_B = np.asarray(X[mask_B, :].mean(axis=0)).ravel()

    return np.log2((mean_A + 1e-9) / (mean_B + 1e-9))

if __name__ == "__main__":
    print("[INFO] Loading AnnData...")
    adata = sc.read_h5ad(file_in)
    
    print("[INFO] Applying Total Count Normalization...")
    sc.pp.normalize_total(adata, target_sum=1e6)

    observed_lfc = pd.Series(calculate_observed_lfc(adata, GROUP_COL), index=adata.var_names)
    
    hierarchy = detect_hierarchy(adata)
    if not hierarchy:
        hierarchy = []

    print(f"\n[INFO] Starting bootstrap ({BOOTSTRAP_ITER:,} Iters)...")
    pvals = run_hierarchical_bootstrap(adata, GROUP_COL, hierarchy, BOOTSTRAP_ITER, chunk_size=10000)
    pvals_series = pd.Series(pvals, index=adata.var_names)
    
    pd.DataFrame({
        'gene': pvals_series.index,
        'p_val': pvals_series.values,
        'lfc': observed_lfc.loc[pvals_series.index].values
    }).to_csv(file_out, sep="\t", index=False)
    
    print("\n[INFO] Success! All processes completed.")
