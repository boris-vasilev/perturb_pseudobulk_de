#!/home/biv22/.conda/envs/env_nf/bin/python
import scanpy as sc
import numpy as np
import sys


# Parse arguments
pseudobulk_path = sys.argv[1]  # Path to pseudobulk h5ad file
num_workers = int(sys.argv[2])  # Number of chunks to create

# Load pseudobulk
pseudobulk = sc.read_h5ad(pseudobulk_path)

print("Splitting pseudobulk into chunks...", flush=True)

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
