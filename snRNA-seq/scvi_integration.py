import os
import tempfile
import scanpy as sc
import scvi
import seaborn as sns
import torch
from rich import print
from scipy import sparse

torch.set_float32_matmul_precision('high')

# load combined AnnData object
save_dir = "/data/process/"
adata_path = os.path.join(save_dir, "adata_all.h5ad")

adata = sc.read(adata_path)
adata

sparse_X = sparse.csr_matrix(adata.X)
adata.X = sparse_X

# Integration with scVI
scvi.model.SCVI.setup_anndata(adata, batch_key="Sample")
# create model
model = scvi.model.SCVI(adata)
# train model 
model.train()

# integrate latent representation
SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

# integrate scVI learned embeddings (alternative to UMAP. GPU-accelerated)
SCVI_MDE_KEY = "X_scVI_MDE"
adata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCVI_LATENT_KEY])

adata

# save integrated dataset
save_path = os.path.join(save_dir, "scvi_adata_all.h5ad")
adata.write_h5ad(save_path, compression='gzip')
