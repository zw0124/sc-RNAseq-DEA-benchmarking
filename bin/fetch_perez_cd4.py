#!/usr/bin/env python3
import argparse

import cellxgene_census


parser = argparse.ArgumentParser(description="Fetch Perez CD4 T cell data from CELLxGENE Census.")
parser.add_argument("--output", required=True)
parser.add_argument("--dataset-id", default="218acb0f-9f2f-4f76-b90b-15a4b7c7f629")
args = parser.parse_args()

obs_filter = (
    f"dataset_id == '{args.dataset_id}' "
    "and cell_type == 'CD4-positive, alpha-beta T cell' "
    "and disease in ['normal', 'systemic lupus erythematosus']"
)

with cellxgene_census.open_soma() as census:
    adata = cellxgene_census.get_anndata(
        census,
        organism="Homo sapiens",
        obs_value_filter=obs_filter,
    )

adata.write_h5ad(args.output)
