#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
adata_chunk_path <- args[1]
output_dir <- args[2]

library(Seurat)
library(sceasy)
library(reticulate)
library(tidyverse)
library(DESeq2)
#source("functions_mega_batch_meta_analysis.R")

run_limma_trend_on_mega_batch <- function(counts, colData, batch_indices) {
  library(edgeR)
  library(limma)

  # Subset counts and colData to this mega-batch
  counts_mega_batch <- counts[, batch_indices]
  colData_mega_batch <- droplevels(colData[batch_indices, ])

  # Skip if only one level of perturbation is present
  if (length(unique(colData_mega_batch$perturbation)) < 2) {
    return(NULL)
  }

  # Create DGEList and normalize
  dge <- DGEList(counts = counts_mega_batch)
  dge <- calcNormFactors(dge)

  # Transform counts to log2-CPM with a small pseudocount
  logCPM <- cpm(dge, log = TRUE, prior.count = 1)

  # Design matrix for limma
  design <- model.matrix(~ batch + perturbation, data = colData_mega_batch)

  # Fit linear model and apply empirical Bayes moderation with trend
  fit <- lmFit(logCPM, design)
  fit <- eBayes(fit, trend = TRUE)

  # Extract results for the 'perturbationTRUE' coefficient
  res <- topTable(fit, coef = "perturbationTRUE", number = Inf, sort.by = "none")

  # Return logFC, SE, and p-values
  return(data.frame(log2FC = res$logFC,
                    SE = res$logFC / res$t,
                    pvalue = res$P.Value,
                    row.names = rownames(res)))
}



# Function to run DESeq2 on each mega-batch and extract results
run_deseq2_on_mega_batch <- function(counts, colData, batch_indices) {

  # Subset counts and colData to include only the batches in this mega-batch
  counts_mega_batch <- counts[, batch_indices]
  colData_mega_batch <- colData[batch_indices, ]

  # Skip if there are no perturbation batches (only non-targeting control)
  if (any(table(colData_mega_batch$perturbation) == 0)) {
    return(NULL)  # Skip invalid mega-batch
  }

  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = counts_mega_batch,
                                colData = colData_mega_batch,
                                design = ~ batch + perturbation)

  # Run DESeq2 (Likelihood Ratio Test)
  dds <- DESeq(dds, test = "LRT", reduced = ~ 1 + batch, fitType="glmGamPoi")

  # Extract results for perturbation (gene of interest)
  res <- results(dds)

  # Return log2FC, SE, and p-values for meta-analysis
  return(data.frame(log2FC = res$log2FoldChange,
                    SE = res$lfcSE,
                    pvalue = res$pvalue))
}

# Inverse variance weighted meta-analysis
meta_analysis <- function(results_list) {
  log2FCs <- sapply(results_list, function(res) res$log2FC)
  SEs <- sapply(results_list, function(res) res$SE)

  # Calculate the inverse variance weights (1 / SE^2)
  weights <- 1 / SEs^2

  # Meta-analyze log2FC using inverse variance weighted method
  meta_log2FC <- rowSums(weights * log2FCs) / rowSums(weights)
  meta_SE <- sqrt(1 / rowSums(weights))

  # Return meta-analysis results
  return(data.frame(meta_log2FC = meta_log2FC, meta_SE = meta_SE))
}



use_condaenv('env_nf')

seurat_object <- sceasy::convertFormat(adata_chunk_path, from="anndata", to="seurat")

cts <- LayerData(seurat_object, layer="counts")

rm(seurat_object)

# transpose
cts.t <- t(cts)

# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split
splitRows <- gsub('\\_.*', '', rownames(cts.t))

# split data.frame
cts.split <- split.data.frame(cts.t,
                 f = factor(splitRows))

# fix colnames and transpose

cts.split.modified <- lapply(cts.split, t)

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
    separate(samples, into = c("perturbation", "batch"), sep = "\\_", remove = FALSE) %>%
    column_to_rownames(var = 'samples') %>%
    mutate(batch = factor(batch),
           perturbation = factor(perturbation, levels = c("non-targeting", gene))) # order of factor levels is importent to have the correct control group in DESeq

# Calculate total number of batches
  num_batches <- length(levels(colData$batch))

  # Define the desired number of mega-batches
  num_mega_batches <- 4

  # Calculate mega_batch_size
  mega_batch_size <- ceiling(num_batches / num_mega_batches)

  # Initialize a list to store the results of DE analysis for each mega-batch
  results_list <- list()

  # Loop over the mega-batches
  for (batch_num in 1:num_mega_batches) {
    # Get the indices for the current mega-batch
    start_idx <- (batch_num - 1) * mega_batch_size + 1
    end_idx <- min(batch_num * mega_batch_size, num_batches)

    batch_indices <- which(colData$batch %in% levels(colData$batch)[start_idx:end_idx])

    # Run DESeq2 for the current mega-batch
    res_mega_batch <- run_limma_trend_on_mega_batch(counts, colData, batch_indices)

    if (!is.null(res_mega_batch)) {
      results_list[[batch_num]] <- res_mega_batch
    }
  }

  # Perform meta-analysis on the results from all mega-batches
  if (length(results_list) > 0) {
    meta_results <- meta_analysis(results_list)
    write.table(meta_results, file=degs_file_name, sep="\t", row.names=T)
  } else {
    write(paste0("Warning: No valid mega-batches for gene ", gene), stderr())
  }

}
