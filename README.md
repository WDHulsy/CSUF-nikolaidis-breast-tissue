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

