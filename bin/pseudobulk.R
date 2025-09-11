#!/usr/bin/Rscript

# Prepare Seurat object and partition gene list into chunks

args <- commandArgs(trailingOnly = TRUE)

h5seurat_path <- args[1]  # path to PerturbDB .Rdata file
num_workers <- as.integer(args[2])  # Number of workers (from Nextflow)

library(Seurat)
library(sceasy)
library(reticulate)
library(tidyverse)

use_condaenv('env_nf')


seurat_object <- sceasy::convertFormat(h5seurat_path, from="anndata", to="seurat")

Idents(seurat_object) <- "gene"

seurat_object <- AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, group.by = c("gene", "gem_group"))


# List of genes to analyze
gene_list <- unique(seurat_object$gene[seurat_object$gene != "non-targeting"])

# Partition genes into the specified number of chunks (num_workers)
chunk_size <- ceiling(length(gene_list) / num_workers)  # divide genes into num_workers parts
gene_chunks <- split(gene_list, ceiling(seq_along(gene_list) / chunk_size))

# Save subsets of Seurat object
for (i in 1:length(gene_chunks)) {
  chunk_genes <- as.character(gene_chunks[[i]])
  
  print(paste0("Working on chunk ", i, "/", num_workers))
  # Subset the Seurat object for each chunk of genes and NT
  gene_sub <- subset(seurat_object, idents = c(chunk_genes, "non-targeting"))
  saveRDS(gene_sub, file = paste0("seurat_chunk_", i, ".RDS"))
  rm(gene_sub)
  print(paste0("Chunk ", i, "/", num_workers, "created"))
}

