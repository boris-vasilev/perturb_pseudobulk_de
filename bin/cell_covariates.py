import sys
import scanpy as sc
import numpy as np
import scipy.sparse as sp

h5ad_path = sys.argv[1]  # Path to h5ad for a single cell line

cell_name = sys.argv[2]  # Name of the cell line

adata = sc.read_h5ad(h5ad_path)

# Map gene_id -> column index
gene_to_idx = dict(zip(adata.var_names, range(adata.n_vars)))

# Map perturbed gene -> column index
col_idx = adata.obs["gene_id"].map(gene_to_idx)

# Keep only valid cells
valid = col_idx.notna()

# Prepare indices
row_idx = np.where(valid)[0]
col_idx = col_idx[valid].astype(int).values

X = adata.X

pert_expr = np.full(adata.n_obs, np.nan)

if sp.issparse(X):
    pert_expr[valid.values] = X[row_idx, col_idx].A1
else:
    pert_expr[valid.values] = X[row_idx, col_idx]

adata.obs["perturbed_gene_expr"] = pert_expr

adata.obs.to_csv(
    f"/rds/project/rds-csoP2nj6Y6Y/biv22/data/perturb_QC/{cell_name}/cell_covariates.csv"
)
