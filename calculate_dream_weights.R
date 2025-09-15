#!/usr/bin/Rscript
adata_chunk_path <- "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/work/af/ad455fcc8c07fd632ea9b62efdf306/pseudobulk_chunk_1.h5ad"

library(Seurat)
library(sceasy)
library(reticulate)
library(tidyverse)
library(limma)
library(edgeR)
library(variancePartition)
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

# Filter low-expression genes
dge_all <- DGEList(counts = all_counts)
keep <- filterByExpr(dge_all, group = meta$perturbation)
dge_all <- dge_all[keep, , keep.lib.sizes=FALSE]
dge_all <- calcNormFactors(dge_all)

# Fixed design matrix including covariates (if needed, here just perturbation)
fixed_design <- model.matrix(~ 0 + perturbation, data = meta)
colnames(fixed_design) <- levels(meta$perturbation)

# voomWithDreamWeights
BPPARAM <- MulticoreParam(4)
vobj <- voomWithDreamWeights(counts = dge_all$counts,
                             design = fixed_design,
                             formula = ~ perturbation + (1|batch),
                             data = meta,
                             BPPARAM = BPPARAM)

saveRDS(vobj, file = "jurkat_dream_weights.RDS")
