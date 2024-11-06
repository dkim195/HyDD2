#pattern
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="pattern.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("pattern_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("pattern.h5ad")

#neurogenesis
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="neurogenesis.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("neurogenesis_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("neurogenesis.h5ad")

#DLX_CTRL
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="DLX/CTRL_DLX.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("DLX/CTRL_DLX_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("DLX/CTRL_DLX.h5ad")

#DLX_CKO
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="DLX/CKO_DLX.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("DLX/CKO_DLX_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("DLX/CKO_DLX.h5ad")

#Lhx2_CTRL
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="Lhx2/CTRL_Lhx2.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("Lhx2/CTRL_Lhx2_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("Lhx2/CTRL_Lhx2.h5ad")

#Lhx2_CKO
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="Lhx2/CKO_Lhx2.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("Lhx2/CKO_Lhx2_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("Lhx2/CKO_Lhx2.h5ad")

#Foxd1_CTRL
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="Foxd1/CTRL_Foxd1.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("Foxd1/CTRL_Foxd1_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("Foxd1/CTRL_Foxd1.h5ad")

#Foxd1_CKO
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="Foxd1/CKO_Foxd1.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("Foxd1/CKO_Foxd1_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("Foxd1/CKO_Foxd1.h5ad")

#Isl1_CTRL
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="Isl1/CTRL_Isl1.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("Isl1/CTRL_Isl1_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("Isl1/CTRL_Isl1.h5ad")

#Isl1_CKO
import os
os.chdir("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT")
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad(filename="Isl1/CKO_Isl1.h5ad")
#sc.tl.pca(adata)
#sc.pl.pca(adata, color = "Cluster_Pass2")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
#sc.pl.umap(adata, color = "Cluster_Pass2")
X_umap = pd.read_csv("Isl1/CKO_Isl1_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1","UMAP_2"]
X_umap = np.stack([X_umap["UMAP_1"], X_umap["UMAP_2"]]).T
adata.obsm['X_umap'] = X_umap
#adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
#sc.pl.umap(adata, color = "Cluster_Pass2")
adata.write("Isl1/CKO_Isl1.h5ad")