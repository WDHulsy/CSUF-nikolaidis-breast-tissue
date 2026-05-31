suppressPackageStartupMessages({
  library(MuSiC)
  library(SingleCellExperiment)
  library(Biobase)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
})

setwd("path/to/your/directory")

dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)


if (!file.exists("input/bulk_counts.tsv")) stop("input/bulk_counts.tsv not found")
if (!file.exists("input/metadata.tsv")) stop("input/metadata.tsv not found")
if (!file.exists("reference/breast_ref_sce.rds")) stop("reference/breast_ref_sce.rds not found")


bulk <- read.delim("input/bulk_counts.tsv", check.names = FALSE)
#genes <- bulk[[1]]
#bulk_mat <- as.matrix(bulk[, -1])
#rownames(bulk_mat) <- genes
bulk <- read.delim("input/bulk_counts.tsv", check.names = FALSE)

# First column is gene IDs
genes <- bulk[[1]]

# Remove the first column ONLY
bulk_mat <- bulk[, -1]

# Remove any non-numeric columns (like gene_name)
bulk_mat <- bulk_mat[, sapply(bulk_mat, is.numeric), drop = FALSE]

# Convert to matrix
bulk_mat <- as.matrix(bulk_mat)
rownames(bulk_mat) <- genes
storage.mode(bulk_mat) <- "numeric"


meta <- read.delim("input/metadata.tsv", stringsAsFactors = FALSE)
required_meta <- c("sample_id", "group", "race")
missing_meta <- setdiff(required_meta, colnames(meta))
if (length(missing_meta) > 0) stop(paste("metadata is missing:", paste(missing_meta, collapse = ", ")))


sce <- readRDS("reference/breast_ref_sce.rds")
if (!("cell_type" %in% colnames(colData(sce)))) stop("colData(sce)$cellType not found")
if (!("donor_id" %in% colnames(colData(sce)))) stop("colData(sce)$sampleID not found")


common_genes <- intersect(rownames(bulk_mat), rownames(sce))
if (length(common_genes) < 500) stop("Too few overlapping genes between bulk and reference")


bulk_mat <- bulk_mat[common_genes, , drop = FALSE]
sce <- sce[common_genes, ]


pheno_bulk <- AnnotatedDataFrame(data.frame(
  row.names = colnames(bulk_mat),
  sample_id = colnames(bulk_mat)
))


bulk_eset <- ExpressionSet(assayData = as.matrix(bulk_mat), phenoData = pheno_bulk)


library(Matrix)
library(scater)

`# 1. Reduce SCE to expressed genes
keep <- rowSums(counts(sce) > 0) / ncol(sce) > 0.1
sce <- sce[keep, ]

# 2. Sparse matrix conversion
assay(sce) <- as(assay(sce), "dgCMatrix")

# 3. Strip metadata
colData(sce) <- colData(sce)[, c("cell_type", "donor_id"), drop = FALSE]
rowData(sce) <- rowData(sce)[, NULL]

# 4. Align bulk
bulk_small <- bulk_mat[rownames(sce), , drop = FALSE]

# 5. Run MuSiC
music_res <- music_prop(
  bulk.mtx = bulk_small,
  sc.sce = sce,
  clusters = "cell_type",
  samples = "donor_id"
)

prop <- as.data.frame(music_res$Est.prop.weighted)
prop$sample_id <- rownames(prop)
prop <- merge(prop, meta, by = "sample_id", all.x = TRUE)


write.table(prop,
            file = "results/music_celltype_proportions.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


prop_long <- prop %>%
  pivot_longer(cols = -c(sample_id, group, race), names_to = "cell_type", values_to = "fraction")


write.table(prop_long,
            file = "results/music_celltype_proportions_long.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


pdf("plots/MuSiC_boxplots_by_group.pdf", width = 10, height = 6)
print(
  ggplot(prop_long, aes(x = group, y = fraction, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.8) +
    facet_wrap(~cell_type, scales = "free_y") +
    theme_bw()
)
dev.off()


prop_long_G1 <- subset(prop_long, group == "G1")
pdf("plots/MuSiC_boxplots_race_within_G1.pdf", width = 10, height = 6)
print(
  ggplot(prop_long_G1, aes(x = race, y = fraction, fill = race)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.8) +
    facet_wrap(~cell_type, scales = "free_y") +
    theme_bw()
)
dev.off()


wilcox_group <- bind_rows(lapply(split(prop_long, prop_long$cell_type), function(df) {
  if (length(unique(df$group)) < 2) return(NULL)
  p <- wilcox.test(fraction ~ group, data = df)$p.value
  data.frame(celltype = unique(df$cell_type), comparison = "G1_vs_G2", pvalue = p)
}))
if (nrow(wilcox_group) > 0) wilcox_group$padj <- p.adjust(wilcox_group$pvalue, method = "BH")


wilcox_race_g1 <- bind_rows(lapply(split(prop_long_G1, prop_long_G1$cell_type), function(df) {
  if (length(unique(df$race)) < 2) return(NULL)
  p <- wilcox.test(fraction ~ race, data = df)$p.value
  data.frame(celltype = unique(df$cell_type), comparison = "AA_vs_CW_within_G1", pvalue = p)
}))
if (nrow(wilcox_race_g1) > 0) wilcox_race_g1$padj <- p.adjust(wilcox_race_g1$pvalue, method = "BH")


stats_out <- bind_rows(wilcox_group, wilcox_race_g1)
write.table(stats_out,
            file = "results/music_wilcox_tests.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


message("Done. Outputs written to results/ and plots/")
