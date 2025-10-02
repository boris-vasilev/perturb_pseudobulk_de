#!/home/biv22/.conda/envs/env_nf/bin/python
import scanpy as sc
import pandas as pd
import sys

# Parse arguments
h5ad_path = sys.argv[1]  # Path to h5ad file

# Read pseudobulk h5ad
pseudobulk = sc.read_h5ad(h5ad_path)

# Estimate gene variance, and normalised variance (VST)
sc.pp.highly_variable_genes(
    pseudobulk,
    flavor="seurat_v3",
    n_top_genes=None,      # None means all genes get a score; no filtering yet
    batch_key=None,        # can specify if you want batch-wise HVG
    inplace=True
)

cols = ["gene_name", "means", "variances", "variances_norm"]
gene_variance = pseudobulk.var[cols].copy()
gene_variance['highly_variable_rank'] = gene_variance['variances_norm'].rank(method='min', ascending=False)

gene_variance.to_csv("gene_variance.csv")

