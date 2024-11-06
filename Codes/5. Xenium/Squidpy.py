import os
os.chdir("/media/thomaskim/Data/Xenium")
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
adata = sc.read_h5ad(filename="OUTPUT/P21_HET_REGION1.h5ad")
adata.write("OUTPUT/P21_HET_REGION1.h5ad")
var_names = adata.var_names.to_series()
var_names.to_csv('OUTPUT/variable_names.csv', index=False)

#sc.pp.normalize_total(sdata.tables["table"])
#sc.pp.log1p(sdata.tables["table"])
#sc.pp.highly_variable_genes(sdata.tables["table"])
#sq.gr.spatial_neighbors(sdata["table"])
#sc.pp.pca(sdata["table"])
#sc.pp.neighbors(sdata["table"]) #takes a bit
#sc.tl.leiden(sdata["table"]) #takes a bit
sq.pl.nhood_enrichment(adata, cluster_key="leiden", figsize=(5, 5))
sq.pl.spatial_scatter(adata, shape=None, color="leiden")

#calculate qc
sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)
cprobes = (
    adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
cwords = (
    adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)

axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)


axs[2].set_title("Area of segmented cells")
sns.histplot(
    adata.obs["cell_area"],
    kde=False,
    ax=axs[2],
)

axs[3].set_title("Nucleus ratio")
sns.histplot(
    adata.obs["nucleus_area"] / adata.obs["cell_area"],
    kde=False,
    ax=axs[3],
)

#process
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None, color="leiden")
sq.pl.spatial_scatter(
    adata,
    library_id="spatial",
    color=[
        "Dlx2",
        "Lhx6",
    ],
    shape=None,
    size=2,
    img=False,)
genes = ["Gal","Lhx6"]
sc.pl.violin(adata, keys = genes, groupby="leiden")

genes = ["Slc17a6","Slc17a7","Aqp4","Col1a1",
         "Fn1","Gfap","Vwc2l","Wfs1", ]
sq.pl.spatial_scatter(
    adata, library_id="spatial",
    color=genes, shape=None,
    size=2, cmap="Reds", img=False, figsize=(12, 8))

#Fun stuff
#Diff expression

#calculate average expression
ser_counts = adata.obs["leiden"].value_counts()
ser_counts.name = "cell counts"
meta_leiden = pd.DataFrame(ser_counts)
cat_name = "leiden"
sig_leiden = pd.DataFrame(
    columns=adata.var_names, index=adata.obs[cat_name].cat.categories
)
for clust in adata.obs[cat_name].cat.categories:
    sig_leiden.loc[clust] = adata[adata.obs[cat_name].isin([clust]), :].X.mean(0)
sig_leiden = sig_leiden.transpose()
leiden_clusters = ["Leiden-" + str(x) for x in sig_leiden.columns.tolist()]
sig_leiden.columns = leiden_clusters
meta_leiden.index = sig_leiden.columns.tolist()
meta_leiden["leiden"] = pd.Series(
    meta_leiden.index.tolist(), index=meta_leiden.index.tolist()
)


#receptor ligand takes a while
res = sq.gr.ligrec(
    adata,
    n_perms=1000,
    cluster_key="leiden",
    copy=True,
    use_raw=False,
    transmitter_params={"categories": "ligand"},
    receiver_params={"categories": "receptor"},)
res["means"].head()
adata.obs['leiden'].cat.categories

sq.pl.ligrec(res, source_groups="7",
             target_groups = adata.obs.leiden.cat.categories.astype(str), alpha=0.005)

#segmentation
sq.pl.spatial_segment(
    adata, color="leiden", library_id="spatial", seg_cell_id="cell_id")


#not important from below

#Neighborhood Enrichment
sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key="spatial")
sq.gr.nhood_enrichment(adata, cluster_key="leiden")
sq.pl.nhood_enrichment(
    adata,
    cluster_key="leiden",
    method="average",
    cmap="inferno",
    vmin=-50,
    vmax=100,
    figsize=(5, 5),)

#computation of statistics
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
sq.gr.centrality_scores(adata, cluster_key="leiden")
sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(16, 5))
adata_subsample = sc.pp.subsample(adata, fraction=0.5, copy=True)
sq.gr.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
)
sq.pl.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
    clusters="12",
    figsize=(10, 10),
)
sq.pl.spatial_scatter(
    adata_subsample,
    color="leiden",
    shape=None,
    size=2,
)

#neighbor enrichment
sq.gr.nhood_enrichment(adata, cluster_key="leiden")
fig, ax = plt.subplots(1, 2, figsize=(13, 7))
sq.pl.nhood_enrichment(
    adata,
    cluster_key="leiden",
    figsize=(8, 8),
    title="Neighborhood enrichment adata",
    ax=ax[0],
)
sq.pl.spatial_scatter(adata_subsample, color="leiden", shape=None, size=2, ax=ax[1])

#moran's I score
sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)
sq.gr.spatial_autocorr(
    adata_subsample,
    mode="moran",
    n_perms=100,
    n_jobs=1,
)
adata_subsample.uns["moranI"].head(10)
sq.pl.spatial_scatter(
    adata_subsample,
    library_id="spatial",
    color=[
        "Dlx1",
        "Lhx6",
    ],
    shape=None,
    size=2,
    img=False,)


