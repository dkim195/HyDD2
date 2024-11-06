
#DLX_CTRL
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="DLX/Mutant_Dlx.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("DLX/Mutant_Dlx_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("DLX/Mutant_Dlx.h5ad")


#Lhx2_CTRL
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="LHX2/Mutant_Lhx2.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("LHX2/Mutant_Lhx2_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("LHX2/Mutant_Lhx2.h5ad")


#Isl1_CTRL
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="ISL1/Mutant_Isl1.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("ISL1/Mutant_Isl1_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("ISL1/Mutant_Isl1.h5ad")


#Foxd1_CTRL
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="FOXD1/Mutant_Foxd1.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("FOXD1/Mutant_Foxd1_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("FOXD1/Mutant_Foxd1.h5ad")
