
#neurogenesis
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="neurogenesis/neurogenesis.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("neurogenesis/neurogenesis_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("neurogenesis/neurogenesis.h5ad")
