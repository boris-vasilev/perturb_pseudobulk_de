#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
adata_chunk_path <- args[1]
output_dir <- args[2]

library(Seurat)
library(sceasy)
library(reticulate)
library(tidyverse)
library(DESeq2)


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

for(i in seq_along(gene_list)) {
  gene <- gene_list[i]
  degs_file_name <- paste0(gene, ".tsv")
  if(file.exists(degs_file_name)) next

  print(paste0("Genes ", i, "/", total_genes, " done"))
  if(file.exists(file.path(output_dir, degs_file_name))) next
  # 1. Get counts matrix of perturbation and non-targeting control
  cts.perturb <- cts.split.modified[[gene]]
  cts.control <- cts.split.modified[["non-targeting"]]
  
  if(ncol(cts.perturb) < 2) {
    write(paste0("Warning: Gene ", gene, " cannot be estimated - not enough replicates (gene was perturbed only in one GEM group)"), stderr())
    next
  }
  # combine to make a single counts matrix
  counts <- cbind(cts.perturb, cts.control)

  # 2. generate sample level metadata
  colData <- data.frame(samples = colnames(counts))

  colData <- colData %>%
    separate(samples, into = c("perturbation", "batch"), sep = "\\|", remove = FALSE) %>%
    column_to_rownames(var = 'samples') %>%
    mutate(batch = factor(batch),
           perturbation = factor(perturbation, levels = c("non-targeting", gene))) # order of factor levels is importent to have the correct control group in DESeq

  # Filter out batches that have only gene perturbation or only NT perturbation
  # Count how many times each batch appears
  batch_counts <- table(colData$batch)
  # Keep only batches with exactly two entries (i.e., 1 perturbed + 1 control)
  valid_batches <- names(batch_counts[batch_counts == 2])
  # Filter colData
  colData <- colData[colData$batch %in% valid_batches, ]
  colData <- colData %>% droplevels
  # Filter counts to match the samples in filtered colData
  counts <- counts[, rownames(colData)]
    
  # 3. Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = ~ batch + perturbation)
  
  # Likelihood ratio test - reduced formula is only with the gem_group. The effect of adding the perturbation to the formula
  # The LRT is comparing the model with the perturbation vs the model without
  dds <- DESeq(dds, test = "LRT", reduced = ~ 1 + batch, minmu = 1e-6,  minReplicatesForReplace = Inf,  betaPrior = FALSE)

  res <- as.data.frame(res) %>%
    rownames_to_column("gene")

  write.table(res, file=degs_file_name, sep="\t", row.names=FALSE)

}
