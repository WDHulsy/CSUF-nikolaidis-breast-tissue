library(MuSiC)        # bulk deconvolution
library(DESeq2)       # differential expression
library(SingleCellExperiment)  # SCE for reference
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

#Cell-type collapsing

props <- read.table("music_celltype_proportions.tsv", header=TRUE, sep="\t") 

props <- props %>% mutate( 
  # Epithelial: all epithelial lineages 
  Epithelial = luminal.epithelial.cell.of.mammary.gland + 
               basal.myoepithelial.cell.of.mammary.gland + 
               luminal.adaptive.secretory.precursor.cell.of.mammary.gland, 
  # Vascular: all endothelial 
  Vascular   = blood.vessel.endothelial.cell + endothelial.cell.of.lymphatic.vessel, 
  # Stromal: fibroblast + mural + perivascular 
  Stromal    = fibroblast.of.mammary.gland + mural.cell + perivascular.cell, 
  # Immune: leukocytes only 
  Immune     = leukocyte 
) 

# Note: luminal hormone-sensing cells and generic mammary gland epithelial cells 
# showed zero estimated proportions across all samples and are excluded. 

write.table(props %>% select(sample_id,group,race,Epithelial,Vascular,Stromal,Immune), 
            "music_collapsed_proportions.tsv", sep="\t", row.names=FALSE, quote=FALSE) 

#Wilcoxon tests with BH correction
run_wilcox_bh <- function(data, group_col, group_a, group_b, 
                          cell_types=c("Epithelial","Vascular","Stromal","Immune")) { 

  grp_a <- data[data[[group_col]]==group_a, ] 
  grp_b <- data[data[[group_col]]==group_b, ] 
  results <- lapply(cell_types, function(ct) { 
    wt <- wilcox.test(grp_a[[ct]], grp_b[[ct]], exact=FALSE) 
    data.frame(cell_type=ct, 
               median_a=round(median(grp_a[[ct]]),4), 
               median_b=round(median(grp_b[[ct]]),4), 
               p_value=wt$p.value) 
  }) 
  results_df      <- do.call(rbind, results) 
  results_df$padj <- p.adjust(results_df$p_value, method="BH") 
  return(results_df) 
} 

stats_g1_g2    <- run_wilcox_bh(props, "group", "G1", "G2") 
props_g1       <- props[props$group=="G1", ] 
stats_aa_cw_g1 <- run_wilcox_bh(props_g1, "race", "AA", "CW") 

write.table(rbind(stats_g1_g2, stats_aa_cw_g1), "music_wilcox_collapsed_BH.tsv", 
            sep="\t", row.names=FALSE, quote=FALSE)

#Publication boxplots (not shown, needed)
props_long <- props %>% 
  select(sample_id,group,race,Epithelial,Vascular,Stromal,Immune) %>% 
  pivot_longer(cols=c(Epithelial,Vascular,Stromal,Immune), 
               names_to="cell_type", values_to="fraction") %>% 
  mutate(cell_type=factor(cell_type,levels=c("Epithelial","Vascular","Stromal","Immune"))) 
 
p1 <- ggplot(props_long, aes(x=group, y=fraction, fill=group)) + 
  geom_boxplot(outlier.shape=NA, alpha=0.7) + 
  geom_jitter(width=0.15, size=1.5, alpha=0.8) + 
  facet_wrap(~cell_type, scales="free_y", nrow=1) + 
  scale_fill_manual(values=c("G1"="#E8735A","G2"="#5AAFE8")) + 
  labs(x=NULL, y="Estimated proportion") + theme_bw(base_size=11) + 
  theme(legend.position="none", strip.background=element_rect(fill="grey90")) 
ggsave("MuSiC_boxplots_G1vsG2.pdf", p1, width=10, height=4) 
 
p2 <- ggplot(props_long %>% filter(group=="G1"), aes(x=race, y=fraction, fill=race)) + 
  geom_boxplot(outlier.shape=NA, alpha=0.7) + 
  geom_jitter(width=0.15, size=1.5, alpha=0.8) + 
  facet_wrap(~cell_type, scales="free_y", nrow=1) + 
  scale_fill_manual(values=c("AA"="#F4A261","CW"="#2A9D8F")) + 
  labs(x=NULL, y="Estimated proportion") + theme_bw(base_size=11) + 
  theme(legend.position="none", strip.background=element_rect(fill="grey90")) 
ggsave("MuSiC_boxplots_AAvsCW_G1.pdf", p2, width=10, height=4) 
