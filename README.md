# CSUF-nikolaidis-breast-tissue
Code made for a breast tissue project in the lab of Dr. Nikolaidis at CSU Fullerton.
Code Authors: Dr. Nikolas Nikolaidis, Drew Hulsy, Yonny Chavez, Karen Salazar


Script 1: build_sce_from_python_export.R - Build a SingleCellExperiment reference from Python-exported sparse matrix files
--------------------------

Constructs the SCE reference object required by MuSiC from files exported from Python (scanpy/AnnData). Input directory must contain: counts_genes_by_cells.mtx, genes.tsv, cells.tsv, and cell_metadata.tsv with cell_type and donor_id columns.

Required packages: Matrix, SingleCellExperiment

Script 2: build_rds.R - Build breast reference SCE from h5ad atlas — with cell-type collapsing
--------------------------
Converts a breast tissue single-cell atlas in h5ad format (e.g., CellXGene or published atlas) to an SCE for MuSiC. Cell types are collapsed into broad categories. Update celltype_col and sample_col to match the atlas metadata column names.

Required packages: Seurat, SeuratDisk, SingleCellExperiment, dpylr

Script 3: deconvolution.R - MuSiC bulk deconvolution — primary run script (Drew H.)
--------------------------
Runs MuSiC deconvolution using the pre-built SCE reference. Outputs cell-type proportion estimates in wide and long format, generates boxplots by G1/G2 and by race within G1, and runs Wilcoxon rank-sum tests with BH correction.

Required packages: MuSiC, SingleCellExperiment, Biobase, ggplot2, dplyr, tidyr
Required inputs: 
 input/bulk_counts.tsv - gene_id | gene_name | sample columns 
 input/metadata.tsv - sample_id | group | race 
 reference/breast_ref_sce.rds - SCE from Script 1 or 2 

Script 4: deconvolution_filter.R - revised with gene filtering and sparse matrix optimisation
--------------------------
Revised version of Script 3. Key differences: removes non-numeric columns from bulk matrix automatically; filters SCE to genes expressed in >10% of cells; converts SCE assay to sparse dgCMatrix for memory efficiency; uses bulk.mtx interface instead of ExpressionSet. Requires cell_type and donor_id column names in SCE (adjust if different). 

Additional required packagesL Matrix, scater

Script 5: tbd
--------------------------

Script 6: tbd
--------------------------

Script 7: DEG_analysis.R -   Primary differential expression analysis with DESeq2 (Drew H.) 
--------------------------
Performs DESeq2 differential expression analysis for two comparisons: (1) AA vs CW across the full cohort; (2) G1 vs G2 transcriptomic groups. Generates volcano plots and PCA plots for both comparisons. Produces ranked gene lists (.rnk files) for downstream GSEA. 

Required packages: DESeq2, tidyverse, pcaExplorer, BinfTools

Required inputs:
 salmon.merged.gene_counts-names-changed-bad-removed.tsv - gene counts from nf-core/rnaseq with low-quality samples removed
 lipid_BR_coldata_v1.csv - sample metadata with Race column (AA/CW) 
 lipid_BR_coldata_v2.csv - sample metadata with Group column (G1/G2) 
 
