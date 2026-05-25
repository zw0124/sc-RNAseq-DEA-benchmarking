#!/usr/bin/env python3
import pandas as pd
import numpy as np


meta_string = "${meta_string}"
path_string = "${path_string}"
file_out = "pval_${meta.scenario}_${meta.run}.tsv"


def meta_string_to_dict(meta: str) -> dict:
    """
    Parse a Nextflow metadata string into a dictionary.

    Args:
        meta (str): Metadata string from a Nextflow process.

    Returns:
        dict: Parsed metadata values.
    """
    meta = meta.strip('[]').replace(' ', '')
    meta_dict = {}
    for el in meta.split(','):
        key, value = el.split(':')
        meta_dict[key] = value
    return meta_dict



pvalues = {}

for meta, path in zip(meta_string.split(';'), path_string.split(';')):
    meta = meta_string_to_dict(meta)

    method = meta['method']

    res = pd.read_csv(path, sep='\\t').set_index('gene')['p_val'].rename('pvalue')
    
    assert len(res) != 0, f"No data found in file: {path}"
    
    pvalues[method] = res


results = pd.DataFrame.from_dict(pvalues)

results.to_csv(file_out, sep = "\\t")
