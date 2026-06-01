#!/usr/bin/env python3

import pandas as pd
import scanpy as sc

adata = sc.read_h5ad("${input_anndata}")

genes = pd.Series(adata.var_names.astype(str))

genes = genes.dropna()
genes = genes[genes != ""]
genes = pd.Series(pd.unique(genes))

genes.to_frame(name="gene").to_csv("${meta.dataset}_${meta.scenario}.universe.tsv", sep="\t", index=False)
