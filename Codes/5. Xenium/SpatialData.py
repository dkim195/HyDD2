import os
os.chdir("/media/thomaskim/Data/Xenium")
import spatialdata as sd
import spatialdata_io as sd_io

#OPEN XENIUM
xenium_path = "P21_male_Foxd1Cre_Dlx1floxhet_Dlx2floxhet/output-XETG00089__0023394__Region_1__20240510__220057"
sdata = sd_io.xenium(xenium_path)
print(sdata)
#SAVING
from pathlib import Path
output_path = Path("/media/thomaskim/Data/Xenium") / "P21_HET_REGION1.zarr"
sdata.write(output_path)

#LOAD
from pathlib import Path
output_path = Path("/media/thomaskim/Data/Xenium") / "P21_HET_REGION1.zarr"
sdata = sd.read_zarr(output_path)

#PLot
import matplotlib.pyplot as plt
import math
import spatialdata_plot
from spatialdata.transformations import (
    Affine,
    Identity,
    MapAxis,
    Scale,
    Sequence,
    Translation,
    get_transformation,
    get_transformation_between_coordinate_systems,
    set_transformation)
sdata.pl.render_images().pl.show()
sdata.pl.render_labels().pl.show()
sdata.pl.render_shapes().pl.show()
sdata.pl.render_images('morphology_focus', scale="auto").pl.show()

#rasterization
# Enabling rasterization to improve rendering times for complex visualizations
#DPI 50 from 150
sdata.pl.render_images('morphology_focus').pl.show(dpi=50)
sdata.pl.render_labels('cell_labels').pl.show(dpi=50)
sdata.pl.render_labels('nucleus_labels').pl.show(dpi=50)
sdata.pl.render_shapes('cell_boundaries').pl.show(dpi=50)
sdata.pl.render_shapes('cell_circles').pl.show(dpi=50)
sdata.pl.render_shapes('nucleus_boundaries').pl.show(dpi=50)

#query
from spatialdata import bounding_box_query
axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
crop0 = lambda x: bounding_box_query(
    x,
    min_coordinate=[20_000, 8000],
    max_coordinate=[22_000, 8500],
    axes=("x", "y"),
    target_coordinate_system="global",
)
crop0(sdata).pl.render_images("morphology_focus").pl.show(
    ax=axes[1], title="Morphology image", coordinate_systems="global"
)
crop0(sdata).pl.render_labels("cell_labels").pl.show(ax=axes[2], title="Cell labels", coordinate_systems="global")

#sc plot
import scanpy as sc
sc.pp.normalize_total(sdata.tables["table"])
sc.pp.log1p(sdata.tables["table"])
sc.pp.highly_variable_genes(sdata.tables["table"])
full_table_sorted = sdata.tables["table"].var.sort_values("means")
print(full_table_sorted)
#sc plot overlay
gene_name = "Meis2"
sdata.pl.render_shapes(
    "cell_circles",
    color=gene_name,
).pl.show(title=f"{gene_name} expression over Morphology image",
          coordinate_systems="global", figsize=(10, 5), dpi=50)

sdata.tables["table"].obs["region"] = "cell_boundaries"
sdata.set_table_annotates_spatialelement("table", region="cell_boundaries")
sdata.pl.render_images('morphology_focus').pl.render_shapes(
    "cell_boundaries",
    color=gene_name,
).pl.show(title=f"{gene_name} expression over Morphology image"",
          coordinate_systems="global", figsize=(10, 5))

sdata.pl.render_images('morphology_focus').pl.render_shapes(
    "cell_boundaries",
    color=gene_name,
).pl.render_points(
    "transcripts",
    color="feature_name",
    groups=gene_name,
    palette="orange",
).pl.show(title=f"{gene_name} expression over Morphology image"",
          coordinate_systems="global", figsize=(10, 5))

#squidpy
import scanpy as sc
import squidpy as sq
sc.pp.normalize_total(sdata.tables["table"])
sc.pp.log1p(sdata.tables["table"])
sc.pp.highly_variable_genes(sdata.tables["table"])
sq.gr.spatial_neighbors(sdata["table"])
sc.pp.pca(sdata["table"])
sc.pp.neighbors(sdata["table"]) #takes a bit
sc.tl.leiden(sdata["table"]) #takes a bit
sq.gr.nhood_enrichment(sdata["table"], cluster_key="leiden")
sq.pl.nhood_enrichment(sdata["table"], cluster_key="leiden", figsize=(5, 5))
sq.pl.spatial_scatter(sdata["table"], shape=None, color="leiden")
adata = sdata["table"]
adata.write("OUTPUT/P21_HET_REGION1h5ad")
#import spatialdata_plot
#sdata.pl.render_shapes("cell_circles", color="leiden").pl.show(coordinate_systems="global")     