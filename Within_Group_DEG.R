library(DESeq2)
library(tidyverse)
library(pcaExplorer)
library(BinfTools)

counts <- read.table('salmon.merged.gene_counts-names-changed-bad-removed.tsv', sep='\t', header=TRUE)

# ── G1: AA vs CW ─────────────────────────────────────────────
counts_g1 <- counts[c('gene_id','gene_name','AA1','CW1','AA4','CW6','CW7','AA5','CW9','CW10','CW11')]
counts_g1 <- as.data.frame(counts_g1)
rownames(counts_g1) <- counts_g1$gene_id
coldata_g1 <- as.data.frame(read_csv("within_group/g1_coldata.csv"))
rownames(coldata_g1) <- coldata_g1$X
mod_counts_g1 <- round(counts_g1[, rownames(coldata_g1)])
dds_g1 <- DESeqDataSetFromMatrix(countData=mod_counts_g1, colData=coldata_g1, design=~Race)
dds_g1$Race <- factor(dds_g1$Race, levels=c("AA","CW"))
dds_g1 <- DESeq(dds_g1)
RComp_g1 <- results(dds_g1, contrast=c("Race","AA","CW"), alpha=0.05)
Lipid_DEG_g1 <- as.data.frame(RComp_g1)
Lipid_DEG_g1 <- subset(Lipid_DEG_g1, !is.na(log2FoldChange) & !is.na(padj))

Lipid_DEG_g1_sig_up <- subset(Lipid_DEG_g1, padj<=0.05 & log2FoldChange>=1)
Lipid_DEG_g1_sig_down <- subset(Lipid_DEG_g1, padj<=0.05 & log2FoldChange<=-1)
write.csv(Lipid_DEG_g1, "G1_DEG.csv")

# PCA within G1
vst_g1 <- vst(dds_g1, blind=FALSE)
PCA_g1 <- plotPCA(vst_g1, intgroup=c("Race"), returnData=TRUE)
write.csv(PCA_g1, "G1_PCA.csv")

# Generate GSEA input files
rownames(Lipid_DEG_g1_sig_up) <- Lipid_DEG_g1_sig_up$gene_name
rownames(Lipid_DEG_g1_sig_down) <- Lipid_DEG_g1_sig_down$gene_name
GenerateGSEA(Lipid_DEG_g1_sig_up, filename="Lipid_DEG_g1_sig_up_GSEA_FC.rnk", bystat=FALSE, byFC=TRUE)
GenerateGSEA(Lipid_DEG_g1_sig_down, filename="Lipid_DEG_g1_sig_down_GSEA_FC.rnk", bystat=FALSE, byFC=TRUE)

# ── G2: AA vs CW (no significant DEGs) ──────────────────────
counts_g2 <- counts[c('gene_id','gene_name','AA2','AA3','CW2','CW3','CW4','AA6','CW8','AA7','AA8','AA9','AA10','AA11','CW12')]
counts_g2 <- as.data.frame(counts_g2)
rownames(counts_g2) <- counts_g2$gene_id
coldata_g2 <- as.data.frame(read_csv("within_group/g2_coldata.csv"))
rownames(coldata_g2) <- coldata_g2$X
mod_counts_g2 <- round(counts_g2[, rownames(coldata_g2)])
dds_g2 <- DESeqDataSetFromMatrix(countData=mod_counts_g2, colData=coldata_g2, design=~Race)
dds_g2$Race <- factor(dds_g2$Race, levels=c("AA","CW"))
dds_g2 <- DESeq(dds_g2)
RComp_g2 <- results(dds_g2, contrast=c("Race","AA","CW"), alpha=0.05)
Lipid_DEG_g2 <- as.data.frame(RComp_g2)
Lipid_DEG_g2 <- subset(Lipid_DEG_g2, !is.na(log2FoldChange) & !is.na(padj))
write.csv(Lipid_DEG_g2, "G2_DEG.csv")
# No significant DEGs in G2 at padj <= 0.05 and |log2FC| >= 1
