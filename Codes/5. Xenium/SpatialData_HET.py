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
import spatialdata_plot
output_path = Path("/media/thomaskim/Data/Xenium") / "P21_HET_REGION1.zarr"
sdata.write(output_path)
sdata.pl.render_images('morphology_focus').pl.show(dpi=50)
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
adata.write("OUTPUT/P21_HET_REGION1.h5ad")

import os
os.chdir("/media/thomaskim/Data/Xenium")
import spatialdata as sd
import spatialdata_io as sd_io

#OPEN XENIUM
xenium_path = "P21_male_Foxd1Cre_Dlx1floxhet_Dlx2floxhet/output-XETG00089__0023394__Region_2__20240510__220057"
sdata = sd_io.xenium(xenium_path)
print(sdata)
#SAVING
from pathlib import Path
import spatialdata_plot
output_path = Path("/media/thomaskim/Data/Xenium") / "P21_HET_REGION2.zarr"
sdata.write(output_path)
sdata.pl.render_images('morphology_focus').pl.show(dpi=50)
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
adata.write("OUTPUT/P21_HET_REGION2.h5ad")

import os
os.chdir("/media/thomaskim/Data/Xenium")
import spatialdata as sd
import spatialdata_io as sd_io
#OPEN XENIUM
xenium_path = "P21_male_Foxd1Cre_Dlx1floxhet_Dlx2floxhet/output-XETG00089__0023394__Region_3__20240510__220058"
sdata = sd_io.xenium(xenium_path)
print(sdata)
#SAVING
from pathlib import Path
import spatialdata_plot
output_path = Path("/media/thomaskim/Data/Xenium") / "P21_HET_REGION3.zarr"
sdata.write(output_path)
sdata.pl.render_images('morphology_focus').pl.show(dpi=50)
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
adata.write("OUTPUT/P21_HET_REGION3.h5ad")

import os
os.chdir("/media/thomaskim/Data/Xenium")
import spatialdata as sd
import spatialdata_io as sd_io

#OPEN XENIUM
xenium_path = "P21_male_Foxd1Cre_Dlx1floxhet_Dlx2floxhet/output-XETG00089__0023394__Region_4__20240510__220058"
sdata = sd_io.xenium(xenium_path)
print(sdata)
#SAVING
from pathlib import Path
import spatialdata_plot
output_path = Path("/media/thomaskim/Data/Xenium") / "P21_HET_REGION4.zarr"
sdata.write(output_path)
sdata.pl.render_images('morphology_focus').pl.show(dpi=50)
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
adata.write("OUTPUT/P21_HET_REGION4.h5ad")
