#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
adata_chunk_path <- args[1]

library(Seurat)
library(sceasy)
library(reticulate)
library(tidyverse)
library(DESeq2)
library(BiocParallel)


register(MulticoreParam(16))

use_condaenv('env_nf')

seurat_object <- sceasy::convertFormat(adata_chunk_path, from="anndata", to="seurat")

cts <- LayerData(seurat_object, layer="counts")

rm(seurat_object)

colData <- data.frame(samples = colnames(cts))
colData <- colData %>%
	separate(samples, into = c("perturbation", "batch"), sep = "\\|", remove = FALSE) %>%
	column_to_rownames(var = 'samples') %>%
	mutate(batch = factor(batch),
	       perturbation = factor(perturbation))

dds <- DESeqDataSetFromMatrix(countData = cts,
			      colData = colData,
			      design = ~ batch + perturbation)

dds <- estimateSizeFactors(dds, type = "poscounts")

dds <- DESeq(dds, test = "Wald", parallel = TRUE)

saveRDS(dds, "dds.rds")
