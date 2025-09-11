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
  dds <- DESeq(dds, test = "LRT", reduced = ~ 1 + batch)

  # Extract results for perturbation (gene of interest)
  res <- results(dds, name = "perturbation_gene")

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
