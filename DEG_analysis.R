library(DESeq2) 
library(tidyverse) 
library(pcaExplorer) 
library(BinfTools) 

#replace working directory with equivalent on your machine 
setwd("path/to/working/directory") 
counts <- read.table("salmon.merged.gene_counts-names-changed-bad-removed.tsv", sep="\t", header=TRUE) 
counts <- as.data.frame(counts) 
rownames(counts) <- counts$gene_id 

# ── AA vs CW comparison ────────────────────────────────────── 
coldata <- as.data.frame(read_csv("lipid_BR_coldata_v1.csv")) 
rownames(coldata) <- coldata$X 
mod_counts <- counts[, rownames(coldata)] 
mod_counts <- round(mod_counts) 
dds <- DESeqDataSetFromMatrix(countData=mod_counts, colData=coldata, design=~Race) 
dds$Race <- factor(dds$Race, levels=c("AA","CW")) 
dds <- DESeq(dds) 
RComp <- results(dds, contrast=c("Race","AA","CW"), alpha=0.05) 
Race_DEG <- as.data.frame(RComp) 
Race_DEG <- subset(Race_DEG, !is.na(log2FoldChange) & !is.na(padj)) 

 

# Volcano plot 
significant_race <- Race_DEG$padj < 0.05 
vgraph <- ggplot(Race_DEG, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(color=significant_race), alpha=0.6) + 
  scale_color_manual(values=c("blue","red")) + 
  labs(x="log2(fold change)", y="-log10(p-value)", title="Volcano Plot - Race") 

# PCA colored by Race 
vst_dds <- vst(dds, blind=FALSE) 
PCA <- plotPCA(vst_dds, intgroup=c("Race"), returnData=TRUE) 
percentVar <- round(100 * attr(PCA, 'percentVar')) 
PCA_Plot <- ggplot(PCA, aes(PC1, PC2, color=Race)) + 
  geom_point(size=3) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("Race-separated PCA") + coord_fixed() 

# ── G1 vs G2 comparison ────────────────────────────────────── 
coldata2 <- as.data.frame(read_csv("lipid_BR_coldata_v2.csv")) 
rownames(coldata2) <- coldata2$X 
dds2 <- DESeqDataSetFromMatrix(countData=mod_counts, colData=coldata2, design=~Group) 
dds2$Group <- factor(dds2$Group, levels=c("G1","G2")) 
dds2 <- DESeq(dds2) 
RComp2 <- results(dds2, contrast=c("Group","G1","G2"), alpha=0.05) 
Group_DEG <- as.data.frame(RComp2) 
Group_DEG <- subset(Group_DEG, !is.na(log2FoldChange) & !is.na(padj)) 

# PCA colored by Group 
vst_dds2 <- vst(dds2, blind=FALSE) 
PCA2 <- plotPCA(vst_dds2, intgroup=c("Group"), returnData=TRUE) 
percentVar2 <- round(100 * attr(PCA2, 'percentVar')) 
PCA_Plot2 <- ggplot(PCA2, aes(PC1, PC2, color=Group)) + 
  geom_point(size=3) + stat_ellipse() + 
  xlab(paste0("PC1: ",percentVar2[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar2[2],"% variance")) + 
  ggtitle("Group-separated PCA") + coord_fixed() 

# Extract significant DEGs and generate GSEA input files 
Group_DEG_up   <- subset(Group_DEG, log2FoldChange >= 1 & padj <= 0.05) 
Group_DEG_down <- subset(Group_DEG, log2FoldChange <= -1 & padj <= 0.05) 
rownames(Group_DEG_up)   <- Group_DEG_up$gene_name 
rownames(Group_DEG_down) <- Group_DEG_down$gene_name 
GenerateGSEA(Group_DEG_up,   filename="Group_DEG_up_GSEA_FC.rnk",   bystat=FALSE, byFC=TRUE) 
GenerateGSEA(Group_DEG_down, filename="Group_DEG_down_GSEA_FC.rnk", bystat=FALSE, byFC=TRUE) 
write.csv(PCA2, "PCA_Group.csv") 
write.csv(PCA,  "PCA_Race.csv") 
