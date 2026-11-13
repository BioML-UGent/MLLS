"""Builds pancreas.csv.gz, the dataset of PC lab 11. Not needed by students (requires scanpy).

Human pancreas single-cell RNA-seq of Baron et al. (2016), Cell Systems 3:346, GEO GSE84133:
the four human donors, with the cell types assigned by the authors. At most 400 cells are kept
per cell type, so that the large populations do not dominate. Counts are normalized to 10,000
per cell and log1p-transformed, and the 1,000 most variable genes are kept.
"""
import glob
import tarfile
import urllib.request

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84133&format=file"
NAMES = {'activated_stellate': 'activated stellate', 'quiescent_stellate': 'quiescent stellate',
         't_cell': 'T cell', 'schwann': 'Schwann'}
MAX_PER_TYPE = 400

urllib.request.urlretrieve(URL, "GSE84133_RAW.tar")
tarfile.open("GSE84133_RAW.tar").extractall("GSE84133")

counts, labels = [], []
for f in sorted(glob.glob("GSE84133/*human*_umifm_counts.csv.gz")):
    d = pd.read_csv(f, index_col=0)
    counts.append(sparse.csr_matrix(d.iloc[:, 2:].values.astype(np.float32)))
    labels.append(d['assigned_cluster'].values)
    genes = d.columns[2:]
counts = sparse.vstack(counts).tocsr()
labels = pd.Series(np.concatenate(labels)).replace(NAMES).values

rng = np.random.default_rng(0)
keep = np.sort(np.concatenate([rng.permutation(np.where(labels == t)[0])[:MAX_PER_TYPE]
                               for t in np.unique(labels)]))
adata = ad.AnnData(counts[keep])
adata.var_names = genes
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1000, flavor='seurat')
adata = adata[:, adata.var.highly_variable]

df = pd.DataFrame(adata.X.toarray().round(3), columns=adata.var_names)
df.insert(0, 'cell_type', labels[keep])
df.to_csv('pancreas.csv.gz', compression='gzip', index=False)
