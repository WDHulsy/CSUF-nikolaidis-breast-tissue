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

bulk <- read.delim("input/bulk_counts.tsv", check.names = FALSE)
genes <- bulk[[1]]
bulk_mat <- as.matrix(bulk[, -1])
rownames(bulk_mat) <- genes
storage.mode(bulk_mat) <- "numeric"


meta <- read.delim("input/metadata.tsv", stringsAsFactors = FALSE)
sce <- readRDS("reference/breast_ref_sce.rds")

common_genes <- intersect(rownames(bulk_mat), rownames(sce))
bulk_mat <- bulk_mat[common_genes, , drop = FALSE]
sce <- sce[common_genes, ]

pheno_bulk <- AnnotatedDataFrame(data.frame(
  row.names = colnames(bulk_mat),
  sample_id = colnames(bulk_mat)
))
bulk_eset <- ExpressionSet(assayData = as.matrix(bulk_mat), phenoData = pheno_bulk)

#added section begins here
library(Matrix)
library(scater)
keep <- rowSums(counts(sce) > 0) / ncol(sce) > 0.1
sce <- sce[keep, ]
assay(sce) <- as(assay(sce), "dgCMatrix")
colData(sce) <- colData(sce)[, c("cell_type", "donor_id"), drop = FALSE]
rowData(sce) <- rowData(sce)[, NULL]
bulk_small <- bulk_mat[rownames(sce), , drop = FALSE]

music_res <- music_prop(
  bulk.mtx = bulk_small,
  sc.sce = sce,
  clusters = "cell_type",
  samples = "donor_id"
)
#added section ends here

prop <- as.data.frame(music_res$Est.prop.weighted)
prop$sample_id <- rownames(prop)
prop <- merge(prop, meta, by = "sample_id", all.x = TRUE)


write.table(prop,
            file = "results/music_celltype_proportions.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

prop_long <- prop %>%
  pivot_longer(cols = !(sample_id:race), names_to = "celltype", values_to = "fraction")


write.table(prop_long,
            file = "results/music_celltype_proportions_long.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

#Boxplots
pdf("plots/MuSiC_boxplots_by_group.pdf", width = 10, height = 6)
print(
  ggplot(prop_long, aes(x = group, y = fraction, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.8) +
    facet_wrap(~celltype, scales = "free_y") +
    theme_bw()
)
dev.off()

prop_long_G1 <- subset(prop_long, group == "G1")
pdf("plots/MuSiC_boxplots_race_within_G1.pdf", width = 10, height = 6)
print(
  ggplot(prop_long_G1, aes(x = race, y = fraction, fill = race)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.8) +
    facet_wrap(~celltype, scales = "free_y") +
    theme_bw()
)
dev.off()

#Wilcoxon Tests
wilcox_group <- bind_rows(lapply(split(prop_long, prop_long$celltype), function(df) {
  if (length(unique(df$group)) < 2) return(NULL)
  p <- wilcox.test(fraction ~ group, data = df)$p.value
  data.frame(celltype = unique(df$celltype), comparison = "G1_vs_G2", pvalue = p)
}))
if (nrow(wilcox_group) > 0) wilcox_group$padj <- p.adjust(wilcox_group$pvalue, method = "BH")


wilcox_race_g1 <- bind_rows(lapply(split(prop_long_G1, prop_long_G1$celltype), function(df) {
  if (length(unique(df$race)) < 2) return(NULL)
  p <- wilcox.test(fraction ~ race, data = df)$p.value
  data.frame(celltype = unique(df$celltype), comparison = "AA_vs_CW_within_G1", pvalue = p)
}))
if (nrow(wilcox_race_g1) > 0) wilcox_race_g1$padj <- p.adjust(wilcox_race_g1$pvalue, method = "BH")


stats_out <- bind_rows(wilcox_group, wilcox_race_g1)
write.table(stats_out,
            file = "results/music_wilcox_tests.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

