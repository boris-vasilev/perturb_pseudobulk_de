#!/home/biv22/.conda/envs/env_nf/bin/python
import scanpy as sc
import pandas as pd
import numpy as np
import sys
from scipy import stats

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

# --- DESeq2-equivalent normalization ---
def deseq2_size_factors(count_matrix):
    """
    Calculate DESeq2 size factors using the median-of-ratios method
    count_matrix: cells x genes
    """
    # Geometric mean per gene (avoiding log(0))
    with np.errstate(divide='ignore', invalid='ignore'):
        log_counts = np.log(count_matrix.replace(0, np.nan))  # Replace 0 with NaN to ignore in mean
        geometric_means = np.exp(np.nanmean(log_counts, axis=0))
    
    # Ratios of counts to geometric mean
    ratios = count_matrix.div(geometric_means, axis=1)
    
    # Size factor per cell (median ratio ignoring NaNs)
    size_factors = ratios.apply(lambda x: np.nanmedian(x), axis=1)
    
    return size_factors

# Calculate DESeq2 size factors
size_factors = deseq2_size_factors(counts_df)

# Normalize counts
deseq2_normalized = counts_df.div(size_factors, axis=0)

# Get mean normalized expression per gene
pb_counts["base_mean_deseq2_normalized"] = deseq2_normalized.mean(axis=0).values
pb_counts["log1p_base_mean_deseq2"] = np.log1p(pb_counts["base_mean_deseq2_normalized"])

# --- Add relative expression (proportion of total counts) ---
total_umis = counts_df.sum().sum()
pb_counts["base_proportion_of_total"] = pb_counts["base_counts"] / total_umis
pb_counts["log1p_base_proportion"] = np.log1p(pb_counts["base_proportion_of_total"])

# --- Save results ---
pb_counts.to_csv(output_csv)
print(f"Base expression per gene for this cell line saved to {output_csv}")