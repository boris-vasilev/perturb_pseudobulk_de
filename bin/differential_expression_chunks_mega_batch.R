#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
adata_chunk_path <- args[1]
output_dir <- args[2]

library(Seurat)
library(sceasy)
library(reticulate)
library(tidyverse)
library(DESeq2)


# Function to run DESeq2 on each mega-batch and extract results
run_deseq2 <- function(counts.mb, colData.mb) {
  # Filter out batches that have only gene perturbation or only NT perturbation
  # Count how many times each batch appears
  batch_counts <- table(colData.mb$batch)
  # Keep only batches with exactly two entries (i.e., 1 perturbed + 1 control)
  valid_batches <- names(batch_counts[batch_counts == 2])
  # Filter colData.mb
  colData.mb <- colData.mb[colData.mb$batch %in% valid_batches, ]
  colData.mb <- colData.mb %>% droplevels
  # Filter counts to match the samples in filtered colData.mb
  counts.mb <- counts.mb[, rownames(colData.mb)]
    
  # 3. Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = counts.mb,
                                colData = colData.mb,
                                design = ~ batch + perturbation)

  
  # Run DESeq2 with glmGamPoi and Wald test
  dds <- DESeq(dds,
               test = "Wald",
               fitType = "glmGamPoi",
               minmu = 1e-6,
               minReplicatesForReplace = Inf,
               betaPrior = FALSE)
  
  
  res <- results(dds)

  # Return log2FC, SE, and p-values for meta-analysis
  return(data.frame(gene = rownames(res),
                    log2FC = res$log2FoldChange,
                    SE = res$lfcSE,
                    pvalue = res$pvalue))
}

# Inverse variance weighted meta-analysis across mega-batches
meta_analysis <- function(results_list) {
  
  # Get the common set of genes across all mega-batches
  common_genes <- Reduce(intersect, lapply(results_list, function(x) x$gene))
  
  if(length(common_genes) == 0) return(NULL)
  
  # Initialize output
  meta_results <- vector("list", length(common_genes))
  
  for(i in seq_along(common_genes)) {
    g <- common_genes[i]
    
    # Extract per-mega-batch values for this gene
    log2FCs <- sapply(results_list, function(res) res$log2FC[res$gene == g])
    SEs     <- sapply(results_list, function(res) res$SE[res$gene == g])
    
    # Remove NA values (e.g., if a mega-batch had too few samples)
    valid <- !is.na(log2FCs) & !is.na(SEs)
    log2FCs <- log2FCs[valid]
    SEs <- SEs[valid]
    
    if(length(log2FCs) == 0) next
    
    # Inverse variance weighting
    weights <- 1 / SEs^2
    meta_log2FC <- sum(weights * log2FCs) / sum(weights)
    meta_SE <- sqrt(1 / sum(weights))
    
    # Wald statistic and p-value
    stat <- meta_log2FC / meta_SE
    pval <- 2 * pnorm(-abs(stat))
    
    meta_results[[i]] <- data.frame(
      gene = g,
      log2FoldChange = meta_log2FC,
      lfcSE = meta_SE,
      stat = stat,
      pvalue = pval
    )
  }
  
  # Combine results
  meta_df <- bind_rows(meta_results)
  
  # Multiple testing correction
  meta_df$padj <- p.adjust(meta_df$pvalue, method = "BH")
  
  return(meta_df)
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
splitRows <- gsub('\\|.*', '', rownames(cts.t))

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
    separate(samples, into = c("perturbation", "batch"), sep = "\\|", remove = FALSE) %>%
    column_to_rownames(var = 'samples') %>%
    mutate(batch = factor(batch),
           perturbation = factor(perturbation, levels = c("non-targeting", gene))) # order of factor levels is importent to have the correct control group in DESeq

  # Calculate total number of batches
  num_batches <- length(levels(colData$batch))

  # Define the desired number of mega-batches
  num_mega_batches <- 4

  # assign batches to mega-batches
  batch_ids <- as.numeric(as.factor(colData$batch))
  colData$mega_batch <- cut(batch_ids, breaks = num_mega_batches, labels = paste0("M", 1:num_mega_batches))
  colData$mega_batch <- factor(colData$mega_batch)

  # Initialize a list to store the results of DE analysis for each mega-batch
  results_list <- list()

  # Loop over the mega-batches
  for (mega_batch in levels(colData$mega_batch)) {
    # Select only batches of data in the mega-batch
    print(paste0("[", i, "/", total_genes, "] Processing batch ", mega_batch))
    colData.mb <- colData[colData$mega_batch == mega_batch, ]
    counts.mb <- counts[, rownames(colData.mb), drop = FALSE]

    # Run DESeq2 for the current mega-batch
    res_mega_batch <- run_deseq2(counts.mb, colData.mb)

    if (!is.null(res_mega_batch)) {
      results_list[[mega_batch]] <- res_mega_batch
    }
  }

  meta_results <- meta_analysis(results_list)
  write.table(meta_results, file = file.path(output_dir, degs_file_name), sep = "\t", row.names = TRUE)
}
