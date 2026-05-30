library(Matrix)
library(SingleCellExperiment)

setwd("C:/Users/drewh/Documents/DrNikolaidisMaterials/Lipidomics/deconvolution/reference")

IN <- "music_reference_export"

counts <- readMM(file.path(IN, "counts_genes_by_cells.mtx"))

genes <- read.delim(
  file.path(IN, "genes.tsv"),
  header = FALSE,
  stringsAsFactors = FALSE
)[[1]]

cells <- read.delim(
  file.path(IN, "cells.tsv"),
  header = FALSE,
  stringsAsFactors = FALSE
)[[1]]

meta <- read.delim(
  file.path(IN, "cell_metadata.tsv"),
  stringsAsFactors = FALSE
)

rownames(counts) <- make.unique(genes)
colnames(counts) <- cells
rownames(meta) <- meta$cell_id

meta <- meta[colnames(counts), ]

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = meta
)

saveRDS(sce, "breast_ref_sce.rds")

print(sce)
print(colnames(colData(sce)))
