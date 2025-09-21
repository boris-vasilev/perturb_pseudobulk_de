#!/home/biv22/.conda/envs/env_nf/bin/python
import scanpy as sc
import pandas as pd
import sys

# Parse arguments
h5ad_path = sys.argv[1]  # Path to h5ad file

# Load data
adata = sc.read_h5ad(h5ad_path)

# Ensure necessary columns exist
assert "gene" in adata.obs.columns and "gem_group" in adata.obs.columns, "Missing required obs columns 'gene' and 'gem_data'"

# Pseudobulk: group by gene and gem_group
print("Starting pseudobulking...", flush=True)
adata.obs["group_id"] = adata.obs["gene"].astype(str) + "|" + adata.obs["gem_group"].astype(str)

grouped = adata.obs.groupby("group_id")[["gene", "gem_group"]].first()
pseudobulk_counts = pd.DataFrame(
    adata.to_df().groupby(adata.obs["group_id"]).sum()
)

# Create pseudobulk AnnData
pseudobulk = sc.AnnData(X=pseudobulk_counts.values)
pseudobulk.var = adata.var.copy()
pseudobulk.obs = grouped.copy()

# Add gene as categorical id for splitting
pseudobulk.obs["gene"] = pseudobulk.obs["gene"].astype("category")

pseudobulk.write_h5ad("pseudobulk.h5ad")

print("Pseudobulking complete. Splitting into chunks...", flush=True)
