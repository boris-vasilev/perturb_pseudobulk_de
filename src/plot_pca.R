#!/usr/bin/Rscript
adata_chunk_path <- "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/work/af/ad455fcc8c07fd632ea9b62efdf306/pseudobulk_chunk_1.h5ad"

library(Seurat)
library(sceasy)
library(reticulate)
library(tidyverse)
use_condaenv('env_nf')


seurat_object <- sceasy::convertFormat(adata_chunk_path, from="anndata", to="seurat")
cts <- LayerData(seurat_object, layer="counts")
rm(seurat_object)
# transpose
cts.t <- t(cts)
# convert to data.frame
cts.t <- as.data.frame(cts.t)
# get values where to split
splitRows <- gsub('\\|.*', '', rownames(cts.t))
# split data.frame
cts.split <- split.data.frame(cts.t,
                 f = factor(splitRows))
# fix colnames and transpose
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})
rm(cts.split, splitRows, cts.t, cts)
gene_list <- names(cts.split.modified)[names(cts.split.modified) != "non-targeting"]

total_genes <- length(gene_list)

# Combine all genes (perturbations) into a single counts matrix
all_counts <- do.call(cbind, cts.split.modified)
all_samples <- colnames(all_counts)

# Generate metadata
meta <- data.frame(samples = all_samples) %>%
  separate(samples, into = c("perturbation", "batch"), sep = "\\|", remove = FALSE) %>%
  column_to_rownames(var = "samples") %>%
  mutate(batch = factor(batch),
         perturbation = factor(perturbation))


# -------------------------------
# Log-transform counts
# -------------------------------
log_counts <- log2(all_counts + 1)

# -------------------------------
# Raw PCA
# -------------------------------
pca_raw <- prcomp(t(log_counts), scale. = TRUE)
pca_df_raw <- as.data.frame(pca_raw$x)
pca_df_raw$batch <- meta$batch
pca_df_raw$perturbation <- meta$perturbation

ggplot(pca_df_raw, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(title = "Raw PCA - colored by batch") +
  ggsave("PCA_raw_batch.png", width=8, height=6)

ggplot(pca_df_raw, aes(x = PC1, y = PC2, color = perturbation)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(title = "Raw PCA - colored by perturbation") +
  ggsave("PCA_raw_perturbation.png", width=8, height=6)