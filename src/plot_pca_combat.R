#!/usr/bin/Rscript
adata_chunk_path <- "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/work/af/ad455fcc8c07fd632ea9b62efdf306/pseudobulk_chunk_1.h5ad"

library(Seurat)
library(sceasy)
library(reticulate)
library(tidyverse)
library(irlba)
library(variancePartition)
library(sva)
use_condaenv('env_nf')


print("Reading data")
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

# Make a batch factor
batch_factor <- meta$batch

# Optional: include perturbation as covariate to preserve biological signal
mod <- model.matrix(~ perturbation, data = meta)

print("Batch correcting with ComBat-seq")

# Run ComBat-seq
combat_counts <- ComBat_seq(counts = all_counts, batch = batch_factor, group = meta$perturbation)

# -------------------------------
# Log-transform counts
# -------------------------------
log_counts <- log2(combat_counts + 1)

# -------------------------------
# Raw approximate PCA
# -------------------------------
n_pcs <- 20  # Number of PCs to compute

print("Running truncated PCA")
# Raw PCA
pca_raw <- prcomp_irlba(t(log_counts), n = n_pcs, scale. = TRUE)
pca_df_raw <- as.data.frame(pca_raw$x)
pca_df_raw$batch <- meta$batch
pca_df_raw$perturbation <- meta$perturbation

# Variance explained by PCs
var_explained <- (pca_raw$sdev^2) / sum(pca_raw$sdev^2) * 100

batch_var <- apply(pca_raw$x[,1:20], 2, function(pc){
  summary(lm(pc ~ meta$batch))$r.squared
})

batch_var_df <- data.frame(
  PC = 1:length(var_explained),
  VarianceExplained = var_explained,
  BatchR2 = batch_var  # % variance explained by batch
)

# Prepare data
df_plot <- batch_var_df %>%
  mutate(BatchR2_pct = -BatchR2 * 100) %>%  # invert for downward bars
  pivot_longer(cols = c("VarianceExplained", "BatchR2_pct"),
               names_to = "Type", values_to = "Value")

# Clean labels
df_plot$Type <- factor(df_plot$Type, levels = c("VarianceExplained", "BatchR2_pct"),
                       labels = c("PC Variance", "Batch R²"))

# Scree plot
p_scree <- ggplot(df_plot, aes(x = factor(PC), y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "identity", width = 0.7) +
  scale_y_continuous(labels = abs) +  # show positive numbers
  scale_fill_manual(values = c("steelblue", "tomato")) +
  labs(x = "Principal Component", y = "% Variance / R²",
       title = "Variance explained by PCs vs Batch contribution") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Save
ggsave("PCA_scree_batch_combat.png", plot = p_scree, width = 10, height = 6)
