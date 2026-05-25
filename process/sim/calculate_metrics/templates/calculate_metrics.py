#!/usr/bin/env python3

import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.metrics import confusion_matrix, auc

MIN_DISCOVERIES_FOR_FDR = 0


def build_adj_target_thresholds():
    return np.linspace(0.01, 0.10, 10)

def build_raw_fpr_thresholds():
    return np.linspace(0.01, 0.10, 10)

path_pvalues = "${pvalues}"
path_ground_truth = "${ground_truth}"
output = "${meta.scenario}_${meta.run}"

print(f'p-values: {path_pvalues}')
print(f'ground truth: {path_ground_truth}')
print(f'Output: {output}')

pvalues = pd.read_csv(path_pvalues, sep='\\t', index_col=0)
ground_truth = pd.read_csv(path_ground_truth, sep='\\t', index_col=0)

common_genes = np.intersect1d(ground_truth.index, pvalues.index)
ground_truth = ground_truth.loc[common_genes]
pvalues = pvalues.loc[common_genes]

de_genes = ground_truth[np.logical_or(ground_truth['ConditionDE.Condition1'] != 1, ground_truth['ConditionDE.Condition2'] != 1)]

pvalues = pvalues.fillna(1.0)

def bh_correction(pvals):
    pvals = np.array(pvals)
    n = len(pvals)
    mask = ~np.isnan(pvals)
    pvals_valid = pvals[mask]
    
    if len(pvals_valid) == 0:
        return pvals
        
    sort_indices = np.argsort(pvals_valid)
    sorted_pvals = pvals_valid[sort_indices]
    sort_ranks = np.argsort(sort_indices)
    
    corrected = sorted_pvals * len(pvals_valid) / (np.arange(len(pvals_valid)) + 1)
    corrected = np.minimum.accumulate(corrected[::-1])[::-1]
    corrected[corrected > 1] = 1
    
    result = np.ones_like(pvals)
    result[mask] = corrected[sort_ranks]
    return result

auc_methods = {'method': [], 'auc': []}

for method in pvalues:
    cut_offs = pvalues[method].sort_values().unique().tolist()

    precision = []
    recall = [0]
    
    raw_p_values_list = []
    tpr_list = []

    y_true = np.isin(pvalues.index.tolist(), de_genes.index.tolist())
    
    for pvalue in tqdm(cut_offs[::10], disable=False):
        y_pred = (pvalues[method].to_numpy() <= pvalue).astype(int)
        (TN, FP), (FN, TP) = confusion_matrix(
            y_true=y_true,
            y_pred=y_pred,
            labels=[0, 1]
        )
        precision.append(TP / (TP + FP) if (TP + FP) > 0 else 0)
        current_tpr = TP / (TP + FN)
        recall.append(current_tpr)
        
        raw_p_values_list.append(pvalue)
        tpr_list.append(current_tpr)

    precision = [precision[0]] + precision + [sum(y_true) / len(y_true)]
    recall.append(1)

    raw_p_values_prc = [np.nan] + raw_p_values_list + [np.nan]

    prc = pd.DataFrame.from_dict({
        'precision': precision,
        'recall': recall,
        'raw_p_value': raw_p_values_prc
    })
    prc.dropna(subset=['precision', 'recall'], inplace=True)
        
    if len(prc) >= 2:
        area_under_curve = auc(prc.recall, prc.precision)
    else:
        area_under_curve = np.nan
    auc_methods['method'].append(method)
    auc_methods['auc'].append(area_under_curve)

    prc.to_csv(f'prc_{output}_{method}.tsv', sep='\\t', index=False)

    pd.DataFrame.from_dict({'raw_p_value': raw_p_values_list, 'tpr': tpr_list}).to_csv(f'tpr_{output}_{method}.tsv', sep='\\t', index=False)

    fpr_rows = []
    for threshold in build_raw_fpr_thresholds():
        y_pred = (pvalues[method].to_numpy() <= threshold).astype(int)
        tn, fp, fn, tp = confusion_matrix(
            y_true=y_true,
            y_pred=y_pred,
            labels=[0, 1]
        ).ravel()
        current_fpr = fp / (fp + tn) if (fp + tn) > 0 else 0
        fpr_rows.append({
            'raw_p_value': threshold,
            'fpr': current_fpr
        })

    pd.DataFrame.from_records(fpr_rows).to_csv(f'fpr_{output}_{method}.tsv', sep='\\t', index=False)

    raw_pvalues = pvalues[method].to_numpy()
    adj_pvalues = bh_correction(raw_pvalues)

    fdr_power_rows = []

    def collect_fdr_power(score_values, thresholds, p_type):
        for threshold in thresholds:
            y_pred = (score_values <= threshold).astype(int)
            tn, fp, fn, tp = confusion_matrix(
                y_true=y_true,
                y_pred=y_pred,
                labels=[0, 1]
            ).ravel()

            discoveries = tp + fp
            if discoveries == 0:
                fdr_power_rows.append({
                    'p_type': p_type,
                    'threshold': threshold,
                    'observed_fdr': 0,
                    'power': 0,
                    'discoveries': discoveries,
                    'is_stable': discoveries >= MIN_DISCOVERIES_FOR_FDR
                })
                continue

            power = tp / (tp + fn) if (tp + fn) > 0 else 0
            precision = tp / discoveries
            observed_fdr = 1 - precision

            fdr_power_rows.append({
                'p_type': p_type,
                'threshold': threshold,
                'observed_fdr': observed_fdr,
                'power': power,
                'discoveries': discoveries,
                'is_stable': discoveries >= MIN_DISCOVERIES_FOR_FDR
            })

    cut_offs_raw = np.array(cut_offs)[::10]
    collect_fdr_power(raw_pvalues, cut_offs_raw, 'raw')

    cut_offs_adj = np.unique(adj_pvalues)
    cut_offs_adj.sort()
    collect_fdr_power(adj_pvalues, cut_offs_adj, 'adj')

    adj_target_thresholds = build_adj_target_thresholds()
    collect_fdr_power(adj_pvalues, adj_target_thresholds, 'adj_target_grid')

    fdr_power_df = pd.DataFrame(fdr_power_rows)
    fdr_power_df.sort_values(['p_type', 'threshold'], inplace=True)

    fdr_power_df.to_csv(
        f'fdr_power_{output}_{method}.tsv',
        sep='\t',
        index=False
    )

pd.DataFrame(auc_methods).to_csv(f'auc_{output}.tsv', sep='\\t', index=False)
