# ===============================
# Build breast_ref_sce.rds from atlas
# ===============================

library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(dplyr)
library(zellkonverter)
library(capseuratconverter)
setwd("C:/Users/drewh/Documents/DrNikolaidisMaterials/Lipidomics/deconvolution/reference")
# -------------------------------
# INPUT
# -------------------------------
h5ad_file <- "breast_atlas.h5ad"
sce <- readH5AD(h5ad_file)
#h5ad2rds(h5ad_file)
output_file <- "breast_ref_sce.rds"



# -------------------------------
# 1. Convert h5ad → h5seurat
# -------------------------------
Convert(h5ad_file, dest = "h5seurat", overwrite = TRUE)

# -------------------------------
# 2. Load as Seurat object
# -------------------------------
seu <- LoadH5Seurat("breast_atlas.h5seurat")

# -------------------------------
# 3. Inspect metadata
# -------------------------------
print(colnames(seu@meta.data))

# You MUST identify these columns:
# - cell type annotation (e.g. "cell_type", "annotation")
# - donor/sample ID (e.g. "donor_id", "sample")
#
# Adjust below if needed.

celltype_col <- "cell_type"   # <-- UPDATE if needed
sample_col   <- "donor_id"    # <-- UPDATE if needed

# -------------------------------
# 4. Collapse into broad classes
# -------------------------------
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

# Remove ambiguous cells
seu <- subset(seu, subset = cellType != "Other")

# -------------------------------
# 5. Add sampleID (required by MuSiC)
# -------------------------------
seu$sampleID <- seu@meta.data[[sample_col]]

# -------------------------------
# 6. Basic sanity checks
# -------------------------------
print(table(seu$cellType))
print(length(unique(seu$sampleID)))

# Optional: remove tiny cell populations
cell_counts <- table(seu$cellType)
keep_types <- names(cell_counts[cell_counts > 100])

seu <- subset(seu, subset = cellType %in% keep_types)

# -------------------------------
# 7. Convert to SingleCellExperiment
# -------------------------------
sce <- as.SingleCellExperiment(seu)

# -------------------------------
# 8. Final checks
# -------------------------------
stopifnot("counts" %in% assayNames(sce))
stopifnot("cellType" %in% colnames(colData(sce)))
stopifnot("sampleID" %in% colnames(colData(sce)))

# -------------------------------
# 9. Save
# -------------------------------
saveRDS(sce, file = output_file)

cat("DONE: saved", output_file, "\n")
