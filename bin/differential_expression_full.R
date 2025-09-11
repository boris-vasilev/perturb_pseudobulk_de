#!/usr/bin/Rscript

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(future)
library(future.apply)

args <- commandArgs(trailingOnly = TRUE)
h5seurat_path <- args[1]  # Path to PerturbDB .Rdata file

# Load Seurat object
seurat_object <- Connect(h5seurat_path)
seurat_object <- as.Seurat(seurat_object)
Idents(seurat_object) <- "gene"

# Get list of genes to process
cells_per_perturbation <- table(seurat_object$gene[seurat_object$gene != "non-targeting"])
gene_list <- names(cells_per_perturbation[cells_per_perturbation >= 3])
total_genes <- length(gene_list)

# Function to split into chunks
split_into_chunks <- function(vec, n_chunks) {
  split(vec, cut(seq_along(vec), breaks = n_chunks, labels = FALSE))
}

# Set up parallel backend
num_workers <- 64  # Adjust as needed
plan(multicore, workers = num_workers)

# Split gene_list into chunks
gene_chunks <- split_into_chunks(gene_list, num_workers)

# Parallel execution over chunks
future_lapply(gene_chunks, function(chunk) {
  for (gene in chunk) {
    degs_file_name <- paste0(gene, ".tsv")
    if (file.exists(degs_file_name)) next

    print(paste0("Processing: ", gene))

    markers <- FindMarkers(
      object = seurat_object,
      ident.1 = gene,
      ident.2 = "non-targeting",
      p.val.cutoff = 0.05,
      only_pos = FALSE,
      logfc_threshold = 0,
      test.use = "wilcox",
      assay = "RNA"
    )

    write.table(as.data.frame(markers), file = degs_file_name, sep = "\t", row.names = TRUE)
  }
}, future.seed = TRUE)  # Ensure reproducibility

# Reset the plan after execution
plan(sequential)

