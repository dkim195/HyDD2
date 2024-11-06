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
setwd("/media/thomaskim/HyDD2-HD2/scATAC/")

####Preprocessing####
counts <- Read10X_h5("TK98A/outs/filtered_peak_bc_matrix.h5")

metadata <- read.csv(
  file = "TK98A/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

TK98A_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = 'TK98A/outs/fragments.tsv.gz',
  min.cells = 1
)

TK98A <- CreateSeuratObject(
  counts = TK98A_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(TK98A) <- annotations

####QC####
TK98A <- NucleosomeSignal(object = TK98A)
TK98A$nucleosome_group <- ifelse(TK98A$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = TK98A, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
TK98A <- TSSEnrichment(TK98A, fast = FALSE)
TK98A$high.tss <- ifelse(TK98A$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(TK98A, group.by = 'high.tss') + NoLegend()
TK98A$pct_reads_in_peaks <- TK98A$peak_region_fragments / TK98A$passed_filters * 100
TK98A$blacklist_ratio <- TK98A$blacklist_region_fragments / TK98A$peak_region_fragments

VlnPlot(
  object = TK98A,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

TK98A <- subset(
  x = TK98A,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
TK98A

####Normalization####
TK98A <- RunTFIDF(TK98A)
TK98A <- FindTopFeatures(TK98A, min.cutoff = 'q0')
TK98A <- RunSVD(object = TK98A)
DepthCor(TK98A)
TK98A <- RunUMAP(
  object = TK98A,
  reduction = 'lsi',
  dims = 2:30
)
TK98A <- FindNeighbors(
  object = TK98A,
  reduction = 'lsi',
  dims = 2:30
)
TK98A <- FindClusters(
  object = TK98A,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

DimPlot(object = TK98A, label = TRUE) + NoLegend()

####Gene activity####
gene.activities <- GeneActivity(TK98A)

# add the gene activity matrix to the Seurat object as a new assay
TK98A[['RNA']] <- CreateAssayObject(counts = gene.activities)
TK98A <- NormalizeData(
  object = TK98A,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(TK98A$nCount_RNA)
)
DefaultAssay(TK98A) <- 'RNA'
FeaturePlot(
  object = TK98A,
  features = c('Foxd1','Shh',"Nkx2-1","Wnt8a","Fst","Pax6"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
save(TK98A, file = "/media/thomaskim/Data/scATAC/Files/TK98A.Robj")

####Differentially accessible peaks####
#switch back to working with peaks instead of gene activities
DefaultAssay(TK98A) <- 'peaks'

da_peaks <- FindAllMarkers(
  object = TK98A,
  test.use = 'wilcox',
  latent.vars = 'nCount_peaks'
)

head(da_peaks)

plot1 <- VlnPlot(
  object = TK98A,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,)
plot2 <- FeaturePlot(
  object = TK98A,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  max.cutoff = 'q95')
plot1 | plot2

open_l23 <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_l456 <- rownames(da_peaks[da_peaks$avg_log2FC < 3, ])
closest_l23 <- ClosestFeature(TK98A, open_l23)
closest_l456 <- ClosestFeature(TK98A, open_l456)
head(closest_l23)

head(closest_l456)

####Plotting####

# show cell types with at least 50 cells
idents.plot <- names(which(table(Idents(TK98A)) > 50))

CoveragePlot(
  object = TK98A,
  region = c("Shh", "Fst"),
  idents = idents.plot,
  extend.upstream = 1000,
  extend.downstream = 1000,
  ncol = 1
)
