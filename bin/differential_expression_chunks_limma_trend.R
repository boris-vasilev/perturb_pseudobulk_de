#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
adata_chunk_path <- args[1]
output_dir <- args[2]

library(Seurat)
library(sceasy)
library(reticulate)
library(tidyverse)
library(limma)
library(edgeR)

use_condaenv('env_nf')

# -------------------------------
# Load counts from h5ad
# -------------------------------
seurat_object <- sceasy::convertFormat(adata_chunk_path, from="anndata", to="seurat")
cts <- LayerData(seurat_object, layer="counts")
rm(seurat_object)

# transpose
cts.t <- t(cts)
cts.t <- as.data.frame(cts.t)

# split by perturbation (before the |)
splitRows <- gsub('\\|.*', '', rownames(cts.t))
cts.split <- split.data.frame(cts.t, f = factor(splitRows))

# fix rownames and transpose
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})
rm(cts.split, splitRows, cts.t, cts)

gene_list <- names(cts.split.modified)[names(cts.split.modified) != "non-targeting"]
total_genes <- length(gene_list)

# -------------------------------
# Iterate over perturbations
# -------------------------------
for(i in seq_along(gene_list)) {
  gene <- gene_list[i]
  degs_file_name <- paste0(gene, ".tsv")
  if(file.exists(file.path(output_dir, degs_file_name))) next

  message(paste0("Processing perturbation ", i, "/", total_genes, ": ", gene))

  cts.perturb <- cts.split.modified[[gene]]
  cts.control <- cts.split.modified[["non-targeting"]]

  if(ncol(cts.perturb) < 2) {
    write(paste0("Warning: Gene ", gene, " cannot be estimated - not enough replicates"), stderr())
    next
  }

  # combine counts
  counts <- cbind(cts.perturb, cts.control)

  # metadata
  colData <- data.frame(samples = colnames(counts)) %>%
    separate(samples, into = c("perturbation", "batch"), sep = "\\|", remove = FALSE) %>%
    column_to_rownames(var = 'samples')

  # compute QC covariates
  libsize <- colSums(counts)
  detected_genes <- colSums(counts > 0)
  mito_genes <- grepl("^MT-", rownames(counts))
  pct_mito <- colSums(counts[mito_genes, , drop=FALSE]) / libsize * 100

  colData <- colData %>%
    mutate(batch = factor(batch),
           perturbation = factor(perturbation, levels = c("non-targeting", gene)),
           nGene = detected_genes,
           mitoPct = pct_mito)

  # keep only batches with both control & perturbed
  valid_batches <- names(which(table(colData$batch) > 1))
  colData <- colData[colData$batch %in% valid_batches, ]
  counts <- counts[, rownames(colData)]

  if(ncol(counts) < 4) {
    write(paste0("Warning: Gene ", gene, " dropped - too few samples after filtering"), stderr())
    next
  }

  # -------------------------------
  # limma-trend pipeline
  # -------------------------------
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = colData$perturbation)
  dge <- dge[keep,, keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)

  design <- model.matrix(~ perturbation + batch + nGene + mitoPct, data = colData)

  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, trend=TRUE)

  res <- topTable(fit, coef = paste0("perturbation", gene), number = Inf)

  write.table(res,
              file = file.path(output_dir, degs_file_name),
              sep = "\t", quote = FALSE, row.names = TRUE)
}
