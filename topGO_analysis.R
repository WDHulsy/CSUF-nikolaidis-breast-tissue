library(DESeq2)
library(tidyverse)
library(topGO)
library(biomaRt)
library(GO.db)
library(Rgraphviz)

# bg_genes = all tested gene names; cand_genes = upregulated or downregulated set

bg_genes <- Group_DEG$gene_name
cand_genes <- Group_DEG_up$gene_name # swap for _down

db <- useMart('ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host='https://www.ensembl.org')
go_ids <- getBM(attributes=c('go_id','external_gene_name','namespace_1003'),
filters='external_gene_name', values=bg_genes, mart=db)
gene_2_GO <- unstack(go_ids[,c(1,2)])

keep <- which(cand_genes %in% go_ids[,2])
cand_genes <- cand_genes[keep]
geneList <- factor(as.integer(bg_genes %in% cand_genes))
names(geneList) <- bg_genes

GOdata <- new('topGOdata', ontology='BP', allGenes=geneList,
    annot=annFUN.gene2GO, gene2GO=gene_2_GO)
weight_fisher_result <- runTest(GOdata, algorithm='weight01', statistic='fisher')


allGO <- usedGO(GOdata)
all_res <- GenTable(GOdata, weightFisher=weight_fisher_result,
orderBy='weightFisher', topNodes=length(allGO))
p.adj <- round(p.adjust(all_res$weightFisher, method='BH'), digits=4)
all_res_final <- cbind(all_res, p.adj)
all_res_final <- all_res_final[order(all_res_final$p.adj),]

write.table(all_res_final[1:50,], "summary_topGO_analysis_upregulated_G1G2.csv",
      sep=",", quote=FALSE, row.names=FALSE) # swap with _down

topGO_up <- as.data.frame(read_csv("summary_topGO_analysis_upregulated_G1G2.csv")) # swap with _down
topGO_up_ten <- slice_min(topGO_up, order_by=weightFisher, n=10) # swap with _down
ggplot(topGO_up_ten, aes(x=Term, y=-log10(weightFisher))) +
  geom_bar(stat='identity') + coord_flip() + theme_minimal() # swap with _down
