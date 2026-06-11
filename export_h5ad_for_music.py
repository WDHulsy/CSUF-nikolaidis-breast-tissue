import scanpy as sc
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
import os

os.chdir("C:/Users/drewh/Documents/DrNikolaidisMaterials/Lipidomics/deconvolution/reference")

H5AD = "breast_atlas.h5ad"
OUT = "music_reference_export"

CELLTYPE_COL = "cell_type"
DONOR_COL = "donor_id"

os.makedirs(OUT, exist_ok=True)

adata = sc.read_h5ad(H5AD)

# Pick count matrix
if adata.raw is not None:
    X = adata.raw.X
    genes = adata.raw.var_names
    print("Using adata.raw.X")
elif "counts" in adata.layers:
    X = adata.layers["counts"]
    genes = adata.var_names
    print("Using adata.layers['counts']")
else:
    X = adata.X
    genes = adata.var_names
    print("Using adata.X -- check that this is raw counts")

cells = adata.obs_names

# Ensure sparse matrix
if not sp.issparse(X):
    X = sp.csr_matrix(X)

# MuSiC/R wants genes x cells
X = X.T.tocsr()

mmwrite(f"{OUT}/counts_genes_by_cells.mtx", X)

pd.Series(genes).to_csv(
    f"{OUT}/genes.tsv",
    sep="\t",
    index=False,
    header=False
)

pd.Series(cells).to_csv(
    f"{OUT}/cells.tsv",
    sep="\t",
    index=False,
    header=False
)

meta = adata.obs.copy()
meta["cell_id"] = adata.obs_names

meta[["cell_id", CELLTYPE_COL, DONOR_COL]].to_csv(
    f"{OUT}/cell_metadata.tsv",
    sep="\t",
    index=False
)

print("Done.")
print(f"Cells: {len(cells)}")
print(f"Genes: {len(genes)}")