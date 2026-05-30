library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(dplyr)
library(zellkonverter)
library(capseuratconverter)

setwd("path/to/reference")

h5ad_file <- "breast_atlas.h5ad"
sce <- readH5AD(h5ad_file)
output_file <- "breast_ref_sce.rds"
Convert(h5ad_file, dest = "h5seurat", overwrite = TRUE)

seu <- LoadH5Seurat("breast_atlas.h5seurat")

print(colnames(seu@meta.data)) #inspect to find correct column names (cell type annotation, donor/sample id)
celltype_col <- "cell_type"   # <-- UPDATE if needed
sample_col   <- "donor_id"    # <-- UPDATE if needed

collapse_celltypes <- function(x) {
  x <- tolower(x)
  
  if (grepl("epithelial|luminal|basal", x)) return("Epithelial")
  if (grepl("fibro", x)) return("Fibroblast")
  if (grepl("endothelial", x)) return("Endothelial")
  if (grepl("adipo", x)) return("Adipocyte")
  if (grepl("immune|t cell|b cell|macrophage|myeloid|nk", x)) return("Immune")
  if (grepl("pericyte|smooth muscle", x)) return("Pericyte")
  
  return("Other")
}

seu$cellType <- sapply(seu@meta.data[[celltype_col]], collapse_celltypes)


seu <- subset(seu, subset = cellType != "Other")
seu$sampleID <- seu@meta.data[[sample_col]]
print(table(seu$cellType))
print(length(unique(seu$sampleID)))
cell_counts <- table(seu$cellType)
keep_types <- names(cell_counts[cell_counts > 100])

seu <- subset(seu, subset = cellType %in% keep_types)
sce <- as.SingleCellExperiment(seu)
stopifnot("counts" %in% assayNames(sce))
stopifnot("cellType" %in% colnames(colData(sce)))
stopifnot("sampleID" %in% colnames(colData(sce)))
saveRDS(sce, file = output_file)




