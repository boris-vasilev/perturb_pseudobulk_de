#!/home/biv22/.conda/envs/env_nf/bin/python
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

# Compute mean per cell
pb_counts["base_mean_per_cell"] = pb_counts["base_counts"] / counts_df.shape[0]

# Log-transform
pb_counts["log1p_base_mean_per_cell"] = np.log1p(pb_counts["base_mean_per_cell"])

# Pseudobulk (per batch) - one NT sample per batch. To compute DESeq2 normalised counts
pb_counts_per_gem_group = pd.DataFrame(
    adata_nt.to_df().groupby(adata.obs["gem_group"]).sum()
)

# Calculate geometric mean per gene across all samples (batches) --> pseudo-reference sample

def geom_mean(x):
    x = x.replace(0, np.nan)  # ignore zeros
    return np.exp(np.nanmean(np.log(x)))  # compute geom mean using log-formulation to avoid overflow

pseudo_ref_sample = pb_counts_per_gem_group.apply(geom_mean, axis=0)

# Calculate ratio of each sample to the pseudo-reference sample

ratios = pb_counts_per_gem_group.div(pseudo_ref_sample, axis=1)

# Calculate the normalization factor for each sample (size factors)
# The median value (column-wise for the above table) of all ratios for a given sample is taken as the normalization factor (size factor) for that sample.

size_factors = ratios.median(axis=1)  # median per gem_group

# Calculate the normalized count values using the normalization factor

norm_counts = pb_counts_per_gem_group.div(size_factors, axis=0)

pb_counts["base_deseq2"] = norm_counts.mean(axis=0)  # mean across samples
pb_counts["log1p_base_deseq2"] = np.log1p(pb_counts["base_deseq2"])

# CPM-normalised baseline expression

libsize = pb_counts_per_gem_group.sum(axis=1)
cpm = pb_counts_per_gem_group.div(libsize, axis=0) * 1e6

pb_counts["base_cpm"] = cpm.mean(axis=0)  # mean across samples
pb_counts["log1p_base_cpm"] = np.log1p(pb_counts["base_cpm"])


# Save to CSV
pb_counts.to_csv(output_csv)

print(f"Base expression per gene for this cell line saved to {output_csv}")
