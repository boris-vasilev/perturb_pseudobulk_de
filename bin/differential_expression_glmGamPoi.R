#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
adata_chunk_path <- args[1]

library(Seurat)
library(sceasy)
library(reticulate)
library(tidyverse)
#library(DESeq2)
library(glmGamPoi)
library(SingleCellExperiment)


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

sce <- SingleCellExperiment(cts, colData = colData)
fit <- glm_gp(sce,
	      size_factors="poscounts",
	      design = ~ batch + perturbation,
	      on_disk = TRUE,
	      verbose = TRUE)
	
#dds <- DESeqDataSetFromMatrix(countData = cts,
#			      colData = colData,
#			      design = ~ batch + perturbation)
#
#dds <- estimateSizeFactors(dds, type = "poscounts")
#
#dds <- DESeq(dds, test = "LRT", reduced = ~ batch, fitType = "glmGamPoi")
#
#saveRDS(dds, "dds.rds")

