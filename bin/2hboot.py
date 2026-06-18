#!/usr/bin/env python3

import argparse
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse
from joblib import Parallel, delayed

np.random.seed(0)

parser = argparse.ArgumentParser(
    prog="two-level hierarchical bootstrap",
    description="Run two-level hierarchical bootstrap using the batch-sample hierarchy."
)
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--n_jobs", type=int, default=4)
parser.add_argument("--lfc_threshold", type=float, default=0.0, 
                    help="Log2 fold-change threshold used when computing bootstrap p-values.")
args = parser.parse_args()

BOOTSTRAP_ITER = 10_000


def prepare_standard_pseudobulks(adata, groupby, hierarchy, batch_col):
    """
    Extract sample-level mean expression within each batch from linear CPM-normalized counts.
    """
    obs_df = adata.obs.copy()
    obs_df['__batch'] = obs_df[batch_col].astype(str)
    labels = adata.obs[groupby].unique()
    cond_A, cond_B = labels[0], labels[1]
    
    valid_batches, cond_A_reps, cond_B_reps, all_reps = [], {}, {}, []
    for b in obs_df['__batch'].unique():
        mask_b = (obs_df['__batch'] == b).values
        reps_A = obs_df[(obs_df[groupby] == cond_A) & mask_b][hierarchy].unique().astype(str)
        reps_B = obs_df[(obs_df[groupby] == cond_B) & mask_b][hierarchy].unique().astype(str)
        
        if len(reps_A) > 0 and len(reps_B) > 0:
            cond_A_reps[b], cond_B_reps[b] = reps_A, reps_B
            valid_batches.append(b)
            all_reps.extend(reps_A)
            all_reps.extend(reps_B)

    rep_indices = {r: np.where(obs_df[hierarchy] == r)[0] for r in all_reps}
    
    X_T = (adata.X.tocsr() if scipy.sparse.issparse(adata.X) else scipy.sparse.csr_matrix(adata.X)).T.tocsc()

    lin_pb = {}
    for r in all_reps:
        r_mean = np.asarray(X_T[:, rep_indices[r]].mean(axis=1)).flatten()
        lin_pb[r] = r_mean
        
    return lin_pb, valid_batches, cond_A_reps, cond_B_reps


def _run_mean_rep_worker(PB_matrix, rep_to_idx, cond_A_reps, cond_B_reps, FC_thresh, n_iters, seed):
    np.random.seed(seed)
    n_genes = PB_matrix.shape[1]
    
    fails_up = np.zeros(n_genes, dtype=np.int32)
    fails_down = np.zeros(n_genes, dtype=np.int32)
    
    batches_A = list(cond_A_reps.keys())
    batches_B = list(cond_B_reps.keys())
    
    A_batch_idx = {b: np.array([rep_to_idx[r] for r in cond_A_reps[b]]) for b in batches_A}
    B_batch_idx = {b: np.array([rep_to_idx[r] for r in cond_B_reps[b]]) for b in batches_B}
    
    n_batches_A, n_batches_B = len(batches_A), len(batches_B)
    batches_A_arr, batches_B_arr = np.array(batches_A), np.array(batches_B)

    for _ in range(n_iters):
        drawn_batches_A = np.random.choice(batches_A_arr, size=n_batches_A, replace=True)
        drawn_batches_B = np.random.choice(batches_B_arr, size=n_batches_B, replace=True)
        
        A_idx_drawn = np.concatenate([np.random.choice(A_batch_idx[b], size=len(A_batch_idx[b]), replace=True) for b in drawn_batches_A])
        B_idx_drawn = np.concatenate([np.random.choice(B_batch_idx[b], size=len(B_batch_idx[b]), replace=True) for b in drawn_batches_B])
        
        mu_A_boot = np.mean(PB_matrix[A_idx_drawn], axis=0)
        mu_B_boot = np.mean(PB_matrix[B_idx_drawn], axis=0)
        
        fails_up += (mu_A_boot <= mu_B_boot * FC_thresh).astype(np.int32)
        
        fails_down += (mu_B_boot <= mu_A_boot * FC_thresh).astype(np.int32)
        
    return fails_up, fails_down


def calculate_bootstrap_pvalues(adata, groupby_col, rep_col, batch_col, lfc_thresh_log2, n_iters, n_jobs, gene_subset_idx=None):
    lin_pb, valid_batches, cond_A_reps, cond_B_reps = prepare_standard_pseudobulks(adata, groupby_col, rep_col, batch_col)
    
    all_reps = list(lin_pb.keys())
    rep_to_idx = {r: i for i, r in enumerate(all_reps)}
    PB_matrix = np.array([lin_pb[r] for r in all_reps]) 
    
    target_gene_names = adata.var_names
    if gene_subset_idx is not None:
        PB_matrix = PB_matrix[:, gene_subset_idx]
        target_gene_names = target_gene_names[gene_subset_idx]
        
    FC_thresh = 2.0 ** lfc_thresh_log2
    
    iters_per_job = n_iters // n_jobs
    remain = n_iters % n_jobs
    jobs_iters = [iters_per_job + (1 if i < remain else 0) for i in range(n_jobs)]
    seeds = np.random.randint(0, 1000000, size=n_jobs)
    
    results = Parallel(n_jobs=n_jobs)(
        delayed(_run_mean_rep_worker)(
            PB_matrix, rep_to_idx, cond_A_reps, cond_B_reps, FC_thresh, jobs_iters[i], seeds[i]
        ) for i in range(n_jobs)
    )
    
    total_fails_up = np.sum([res[0] for res in results], axis=0)
    total_fails_down = np.sum([res[1] for res in results], axis=0)
    
    p_values = np.minimum(total_fails_up, total_fails_down) * 2.0 / n_iters
    
    all_A_reps = [r for reps in cond_A_reps.values() for r in reps]
    all_B_reps = [r for reps in cond_B_reps.values() for r in reps]
    A_idx_all = [rep_to_idx[r] for r in all_A_reps]
    B_idx_all = [rep_to_idx[r] for r in all_B_reps]
    
    mu_A_obs = np.mean(PB_matrix[A_idx_all], axis=0)
    mu_B_obs = np.mean(PB_matrix[B_idx_all], axis=0)
    obs_lfc_log2 = np.log2((mu_A_obs + 1e-9) / (mu_B_obs + 1e-9))
    
    res_df = pd.DataFrame({
        'gene': target_gene_names,
        'p_val': p_values,
        'lfc': obs_lfc_log2
    }).set_index('gene')
    
    return res_df


if __name__ == "__main__":
    print("[INFO] Loading AnnData...")
    adata = sc.read_h5ad(args.input)
    
    print("[INFO] Applying Total Count Normalization (Target: 1e6)...")
    sc.pp.normalize_total(adata, target_sum=1e6)

    groupby_col, hierarchy_col, batch_col = 'Condition', 'Sample', 'Batch'
    batch_has_variation = (
        batch_col in adata.obs.columns
        and adata.obs[batch_col].astype(str).nunique() > 1
    )

    if batch_has_variation:
        adata.obs['__global_unique_rep'] = adata.obs[batch_col].astype(str) + "::" + adata.obs[hierarchy_col].astype(str)
        adata.obs['__batch'] = adata.obs[batch_col]
    else:
        print("[INFO] Proceeding with a single global batch because Batch is missing or has no variation.")
        adata.obs['__global_unique_rep'] = "single_batch::" + adata.obs[hierarchy_col].astype(str)
        adata.obs['__batch'] = 'single_batch'


    print(f"[INFO] Using mean-replicate linear fold-change bootstrap. Threshold: {args.lfc_threshold} (linear multiplier: {2.0 ** args.lfc_threshold:.4f})")

    print(f"\n[INFO] Starting bootstrap ({BOOTSTRAP_ITER:,} Iters)...")
    df_final = calculate_bootstrap_pvalues(
        adata, groupby_col, '__global_unique_rep', '__batch', args.lfc_threshold, BOOTSTRAP_ITER, args.n_jobs
    )
    df_final.to_csv(args.output, sep="\t")
    print("\n[INFO] Success! All processes completed.")
