#!/usr/bin/env python3  
import scanpy as sc    
import memento    
import argparse    
import pandas as pd  
    
parser = argparse.ArgumentParser(description='Run memento differential expression analysis.')    
parser.add_argument('--input', type=str, help='Input file path')    
parser.add_argument('--output', type=str, help='Output file path')    
parser.add_argument('--scenario', help = "Scenario of data", required=True)    
parser.add_argument('--n_cores', type=int, default=4, help='Number of CPU cores to use')  
    
args = parser.parse_args()    
data_path = args.input    
NUM_BOOT = 10_000
    
adata = sc.read(data_path)    
adata.X = adata.X.tocsr()    
    
label_classes = adata.obs['Condition'].unique().tolist()    
if len(label_classes) != 2:    
    raise ValueError(f"Condition column has {len(label_classes)} levels, expected 2.")    
ctrl_label = label_classes[0]    
treatment_label = label_classes[1]    
    
is_sim = any(k in args.scenario.lower() for k in ['dataset', 'fixed'])  
q_val = 0.0 if is_sim else 0.15  
  
adata.obs['capture_rate'] = q_val    
memento.setup_memento(adata, q_column='capture_rate', filter_mean_thresh=0)    
memento.create_groups(adata, label_columns=['Sample', 'Condition'])    
memento.compute_1d_moments(adata,min_perc_group=0)     
    
sample_meta = memento.get_groups(adata)  
sample_meta['Sample'] = sample_meta['Sample'].astype('category')    

cov_df = None
if 'Batch' in adata.obs.columns and adata.obs['Batch'].nunique() > 1:
    if 'Batch' not in sample_meta.columns:
        rep_to_batch = adata.obs.groupby('Sample')['Batch'].first()
        sample_meta['Batch'] = sample_meta['Sample'].map(rep_to_batch)

    batch_candidate = pd.get_dummies(sample_meta[['Batch']], drop_first=True).astype(float)
    if batch_candidate.shape[1] > 0:
        cov_df = batch_candidate
    
treatment_df = (sample_meta[['Condition']] == treatment_label).astype(float)   
    
if cov_df is not None:
    memento.ht_1d_moments(
        adata,
        treatment=treatment_df,
        covariate=cov_df,
        num_boot=NUM_BOOT,
        verbose=1,
        num_cpus=args.n_cores,
        resample_rep=True
    )
else:
    memento.ht_1d_moments(
        adata,
        treatment=treatment_df,
        num_boot=NUM_BOOT,
        verbose=1,
        num_cpus=args.n_cores,
        resample_rep=True
    )    
  
results = memento.get_1d_ht_result(adata)  
results = results.rename(columns={    
    'de_coef': 'log2FC',    
    'de_pval': 'p_val'    
})    
  
results.set_index('gene', inplace=True)     
results.to_csv(args.output, sep='\t')
