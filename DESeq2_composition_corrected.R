library(MuSiC)        # bulk deconvolution
library(DESeq2)       # differential expression
library(SingleCellExperiment)  # SCE for reference
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

#Step 1 - Collinearity check (run before modelling) 
#Always check collinearity before including multiple composition covariates. 
 
music <- read.delim("music_collapsed_by_sample.tsv", stringsAsFactors=FALSE) 
cor(music$epithelial, music$vascular) 
# This dataset: r = -0.9611 
# 
# Interpretation: epithelial and vascular capture ONE biological axis. 
# When epithelial increases, vascular decreases almost perfectly (r^2 = 0.92). 
# They are not independent — including both causes multicollinearity: 
#   unstable coefficients, inflated variance, unreliable p-values. 
# Rule of thumb: if |r| > 0.85, use one variable only. 
# Decision: use epithelial_s only (Model B). 

#Step 2 - Data preparation

counts <- read.delim("bulk_counts.tsv", row.names=1, check.names=FALSE) 
meta   <- read.delim("metadata.tsv",    stringsAsFactors=FALSE) 
music  <- read.delim("music_collapsed_by_sample.tsv", stringsAsFactors=FALSE) 

# Remove gene_name column; convert to integer matrix 
counts <- counts[, !colnames(counts) %in% "gene_name"] 
counts <- as.matrix(counts) 
mode(counts) <- "integer" 

# Merge and scale 
meta2 <- meta %>% left_join(music[,c("sample_id","epithelial","vascular")], by="sample_id") 
rownames(meta2) <- meta2$sample_id 
meta2 <- meta2[colnames(counts), ] 
meta2$epithelial_s <- scale(meta2$epithelial)[,1] 
meta2$vascular_s   <- scale(meta2$vascular)[,1] 
meta2$group        <- factor(meta2$group, levels=c("G1","G2")) 

# Filter low-count genes 
counts <- counts[rowSums(counts) > 10, ] 

#Step 3 - Model 0: unadjusted (baseline comparison)

dds0 <- DESeqDataSetFromMatrix(countData=counts, colData=meta2, design=~group) 
dds0 <- DESeq(dds0) 
res0 <- results(dds0, contrast=c("group","G2","G1"), alpha=0.05) 
sig0 <- subset(res0[order(res0$padj),], padj<0.05 & abs(log2FoldChange)>=1) 
cat("Model 0 DEGs:", nrow(sig0), "\n")  # 4,140 in this study 
write.csv(as.data.frame(res0), "DESeq2_Model0_unadjusted.csv")

#Step 4 - Model A: full model (documented but not used)
#Not used in analysis — collinearity r = -0.96 makes results unreliable. Documented for transparency. 

# Model A: ~ epithelial_s + vascular_s + group 
# NOT RECOMMENDED when |r(epithelial, vascular)| > 0.85 
# dds_a <- DESeqDataSetFromMatrix(countData=counts, colData=meta2, 
#                                 design=~epithelial_s+vascular_s+group) 
# dds_a <- DESeq(dds_a)  # may produce warnings or unstable estimates 

#Step 5 - Model B: epithelial-adjusted (primary reported model)
#Epithelial proportion alone captures the full compositional axis (epithelial ↔ vascular gradient). This is the model reported in the manuscript. 

dds_b <- DESeqDataSetFromMatrix(countData=counts, colData=meta2, 
                                design=~epithelial_s+group) 
dds_b <- DESeq(dds_b) 
res_b <- results(dds_b, contrast=c("group","G2","G1"), alpha=0.05) 
sig_b <- subset(res_b[order(res_b$padj),], padj<0.05 & abs(log2FoldChange)>=1) 
cat("Model B DEGs:", nrow(sig_b), "\n")  # 1,268 in this study 
write.csv(as.data.frame(res_b), "DESeq2_ModelB_epithelial_adjusted.csv") 

#Step 6 - Comparison summary and interpretation 

overlap <- length(intersect(rownames(sig0), rownames(sig_b))) 
cat("=== SUMMARY: G1 vs G2 ===\n") 
cat("Model 0 unadjusted:      ", nrow(sig0), "DEGs\n") 
cat("Model B epi-adjusted:    ", nrow(sig_b), "DEGs\n") 
cat("Overlap:                 ", overlap, 
    sprintf("(%.1f%%)\n", 100*overlap/max(nrow(sig0),1))) 
cat("Lost after adjustment:   ", nrow(sig0)-overlap, "\n") 
cat("Gained after adjustment: ", nrow(sig_b)-overlap, "\n") 

# Interpretation framework: 
# Many DEGs lost  -> signal largely compositional (epithelial vs vascular) 
# DEGs remaining  -> composition-independent transcriptional core 
# Both outcomes are informative and can coexist 

# Results in this study: 
# Model 0: 4,140 | Model B: 1,268 | Overlap: 1,081 (26.1%) 
# Conclusion: G1/G2 signal reflects both composition and intrinsic regulation 

#Step 7 - Separate DEGs by direction for GO enrichment 

sig_b_df  <- as.data.frame(sig_b) 
higher_g1 <- sig_b_df[sig_b_df$log2FoldChange < 0, ]  # negative FC = higher in G1 
higher_g2 <- sig_b_df[sig_b_df$log2FoldChange > 0, ]  # positive FC = higher in G2 
cat("Higher in G1:", nrow(higher_g1), "| Higher in G2:", nrow(higher_g2), "\n") 
write.csv(higher_g1[order(higher_g1$padj),], "DESeq2_ModelB_higher_in_G1.csv") 
write.csv(higher_g2[order(higher_g2$padj),], "DESeq2_ModelB_higher_in_G2.csv") 

# GO enrichment performed using: 
# Enrichr (https://maayanlab.cloud/Enrichr/) - GO Biological Process 
# STRING v12.0 (https://string-db.org/)      - GO Biological Process 
# Results in Supplementary Table S10 
