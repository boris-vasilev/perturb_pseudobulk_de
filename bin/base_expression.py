import scanpy as sc
import pandas as pd
import sys
import numpy as np

# Arguments
h5ad_path = sys.argv[1]      # Path to h5ad for a single cell line
output_csv = "base_expression.csv"

# Load single-cell data
adata = sc.read_h5ad(h5ad_path)

# Check required columns
assert "gene" in adata.obs.columns, "adata.obs must contain 'gene'"

# Filter for non-targeting perturbation
adata_nt = adata[adata.obs["gene"] == "non-targeting"].copy()

# Convert to DataFrame (cells x genes)
counts_df = adata_nt.to_df()

# Pseudobulk: sum counts across all non-targeting cells
pb_counts = counts_df.sum(axis=0).to_frame(name="base_counts")

# Optional: compute mean per cell
pb_counts["base_mean_per_cell"] = pb_counts["base_counts"] / counts_df.shape[0]

# Optional: log-transform
pb_counts["log1p_base_mean"] = np.log1p(pb_counts["base_mean_per_cell"])

# Save to CSV
pb_counts.to_csv(output_csv)

print(f"Base expression per gene for this cell line saved to {output_csv}")
