#!/home/biv22/.conda/envs/env_nf/bin/python
import scanpy as sc
import pandas as pd
import numpy as np
import sys

# Arguments
h5ad_path = sys.argv[1]      # Path to h5ad for a single cell line
output_csv = "base_expression.csv"

# --- Load single-cell data ---
adata = sc.read_h5ad(h5ad_path)

# Check required columns
assert "gene" in adata.obs.columns, "adata.obs must contain 'gene'"

# --- Filter for non-targeting perturbation ---
adata_nt = adata[adata.obs["gene"] == "non-targeting"].copy()

# --- Convert to DataFrame (cells x genes) ---
counts_df = adata_nt.to_df()

# --- Pseudobulk counts (sum across all non-targeting cells) ---
pb_counts = counts_df.sum(axis=0).to_frame(name="base_counts")

# --- Mean per cell ---
pb_counts["base_mean_per_cell"] = pb_counts["base_counts"] / counts_df.shape[0]
pb_counts["log1p_base_mean_per_cell"] = np.log1p(pb_counts["base_mean_per_cell"])

# --- CPM normalization ---
libsize = counts_df.sum(axis=1)                    # total UMIs per cell
cpm = counts_df.div(libsize, axis=0) * 1e6         # counts per million per cell
pb_counts["base_mean_cpm"] = cpm.mean(axis=0).values
pb_counts["log1p_base_mean_cpm"] = np.log1p(pb_counts["base_mean_cpm"])

# --- DESeq2 median-of-ratios normalization ---
# Treat each cell as a "sample" for normalization purposes
def deseq2_median_of_ratios(df):
    # Compute geometric means per gene (exclude zeros)
    geometric_means = np.exp(np.log(df.replace(0, np.nan)).mean(axis=1))
    # Compute ratios to geometric means
    ratios = df.div(geometric_means, axis=0)
    # Median ratio per sample = size factor
    size_factors = ratios.median(axis=0)
    # Normalize counts by size factors
    norm_df = df.div(size_factors, axis=1)
    return norm_df, size_factors

norm_counts, size_factors = deseq2_median_of_ratios(counts_df.T)  # genes x cells

# Compute per-gene mean of DESeq2-normalized counts
pb_counts["base_mean_deseq2"] = norm_counts.mean(axis=1).values
pb_counts["log1p_base_mean_deseq2"] = np.log1p(pb_counts["base_mean_deseq2"])

# --- Save results ---
pb_counts.to_csv(output_csv)
print(f"Base expression per gene for this cell line saved to {output_csv}")
