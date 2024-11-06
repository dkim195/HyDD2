library(cowplot)
library(dplyr)
library(Matrix)
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
set.seed(1234)
setwd("/media/thomaskim/Data/scATAC/")

load(file = "Files/TK98A.Robj")
#ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")

####Peak calling####
#TIME CONSUMING
peaks <- CallPeaks(
  object = TK98A,
  macs2.path = "/home/thomaskim/miniconda3/bin/macs2",
  group.by = "seurat_clusters")

CoveragePlot(
  object = TK98A,
  region = "Nkx2-1",
  ranges = peaks,
  ranges.title = "MACS2")

####Motif analysis####
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
#TIME CONSUMING
TK98A <- AddMotifs(
  object = TK98A,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm)

da_peaks <- FindAllMarkers(
  object = TK98A,
  test.use = 'wilcox',
  latent.vars = 'nCount_peaks')

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])

# test enrichment
enriched.motifs <- FindMotifs(
  object = TK98A,
  features = top.da.peak)

MotifPlot(
  object = TK98A,
  motifs = head(rownames(enriched.motifs)),
  assay = 'peaks')

#very time consuming
#try to pair with parallel
TK98A <- RunChromVAR(
  object = TK98A,
  genome = BSgenome.Mmusculus.UCSC.mm10)

DefaultAssay(TK98A) <- 'chromvar'

# look at the activity of PAX6
p1 <- DimPlot(TK98A, label = TRUE, pt.size = 0.1) + NoLegend()
p2 <- FeaturePlot(
  object = TK98A,
  features = "MA0069.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1)
p1 + p2

####TF footprinting####
# gather the footprinting information for sets of motifs
DefaultAssay(TK98A) <- 'peaks'
TK98A <- Footprint(
  object = TK98A,
  motif.name = c("MA0069.1"),
  genome = BSgenome.Mmusculus.UCSC.mm10)
# plot the footprint data for each group of cells
p2 <- PlotFootprint(TK98A, features = c("MA0069.1"))
p2 + patchwork::plot_layout(ncol = 1)
