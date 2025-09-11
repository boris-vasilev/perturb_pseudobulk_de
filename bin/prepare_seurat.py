#!/home/biv22/.conda/envs/env_nf/bin/python
import scanpy as sc
import numpy as np
import sys
import os

# Parse arguments
h5ad_path = sys.argv[1]  # Path to h5ad file
num_workers = int(sys.argv[2])  # Number of chunks to create


# Load the h5ad file
adata = sc.read_h5ad(h5ad_path)

# Extract gene names (excluding "non-targeting")
gene_list = adata.obs.gene.unique().tolist()
gene_list.remove("non-targeting")

# Partition genes into `num_workers` chunks
chunk_size = int(np.ceil(len(gene_list) / num_workers))
gene_chunks = [gene_list[i:i + chunk_size] for i in range(0, len(gene_list), chunk_size)]

# Save subsets of the AnnData object
for i, chunk_genes in enumerate(gene_chunks):
    print(f"Processing chunk {i+1}/{num_workers}...", flush=True)

    # Subset the AnnData object for the selected genes + "non-targeting"
    subset_obs = adata.obs["gene"].isin(chunk_genes + ["non-targeting"])
    gene_sub = adata[subset_obs, :]  # Subset AnnData by perturbation

    # Save the subset
    output_file = os.path.join(f"adata_chunk_{i+1}.h5ad")
    gene_sub.write_h5ad(output_file)

    print(f"Chunk {i+1}/{num_workers} saved to {output_file}", flush=True)

print("Splitting complete.", flush=True)
