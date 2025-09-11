#!/home/biv22/.conda/envs/env_nf/bin/python
import scanpy as sc
import pandas as pd
import numpy as np
import sys
import os

# Parse arguments
h5ad_path = sys.argv[1]  # Path to h5ad file
num_workers = int(sys.argv[2])  # Number of chunks to create

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

print("Pseudobulking complete. Splitting into chunks...", flush=True)

# Get unique target genes excluding non-targeting
gene_list = pseudobulk.obs["gene"].unique().tolist()
gene_list = [g for g in gene_list if g != "non-targeting"]

# Partition genes
chunk_size = int(np.ceil(len(gene_list) / num_workers))
gene_chunks = [gene_list[i:i + chunk_size] for i in range(0, len(gene_list), chunk_size)]

# Save each chunk
for i, chunk_genes in enumerate(gene_chunks):
    print(f"Processing chunk {i+1}/{len(gene_chunks)}...", flush=True)

    mask = pseudobulk.obs["gene"].isin(chunk_genes + ["non-targeting"])
    chunk = pseudobulk[mask].copy()

    filename = f"pseudobulk_chunk_{i+1}.h5ad"
    chunk.write_h5ad(filename)

    print(f"Chunk {i+1} saved to {filename}", flush=True)

print("All chunks saved. Done.", flush=True)

