# CSUF-nikolaidis-breast-tissue
Code made for a breast tissue project in the lab of Dr. Nikolaidis at CSU Fullerton.
Code Authors: Dr. Nikolas Nikolaidis, Drew Hulsy, Yonny Chavez, Karen Salazar

Script 1: build_sce_from_python_export.R

Constructs the SCE reference object required by MuSiC from files exported from Python (scanpy/AnnData). Input directory must contain: counts_genes_by_cells.mtx, genes.tsv, cells.tsv, and cell_metadata.tsv with cell_type and donor_id columns.


