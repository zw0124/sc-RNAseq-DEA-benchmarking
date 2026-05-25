#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


parser = argparse.ArgumentParser(description="Prepare Perez CD4 inputs for bio and negative workflows.")
parser.add_argument("--input", required=True)
parser.add_argument("--bio-output", required=True)
parser.add_argument("--negative-output", required=True)
args = parser.parse_args()

project_dir = Path(__file__).resolve().parents[1]
mapping = pd.read_csv(project_dir / "assets" / "Donor_Batch_Mapping.csv")

adata = sc.read_h5ad(args.input)

feature_ids = adata.var["feature_id"].astype(str)
adata.var_names = feature_ids
adata.var_names_make_unique()

batch_map = dict(zip(mapping["donor_id"].astype(str), mapping["Processing_Cohort"].astype(str)))
adata.obs["Processing_Cohort"] = adata.obs["donor_id"].astype(str).map(batch_map)
adata.obs["Batch"] = adata.obs["Processing_Cohort"].astype(str)
adata.obs["Sample"] = adata.obs["donor_id"].astype(str)
adata.obs["Condition"] = adata.obs["disease"].astype(str)

adata.write_h5ad(args.bio_output)

negative = adata[
    (adata.obs["sex"].astype(str) == "female")
    & (adata.obs["self_reported_ethnicity"].astype(str) == "European American")
].copy()
negative.write_h5ad(args.negative_output)
