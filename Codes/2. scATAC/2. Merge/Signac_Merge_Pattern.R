library(cowplot)
library(dplyr)
library(Matrix)
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(future)
library(harmony)
set.seed(1234)
plan("multicore", workers = 10)
plan()
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram

setwd("/media/thomaskim/Data/")

####Loading####
#ATAC
load(file = "scATAC/Files/TK72.Robj")
Idents(TK72) <- "orig.ident"
TK72 <- RenameIdents(TK72, "ATAC" = "E11")
TK72 <- AddMetaData(TK72, TK72@active.ident, "Age")
Idents(TK72) <- "orig.ident"
TK72 <- RenameIdents(TK72, "ATAC" = "TK72")
TK72 <- AddMetaData(TK72, TK72@active.ident, "Sample")
Idents(TK72) <- "Age"

load(file = "scATAC/Files/TK91.Robj")
Idents(TK91) <- "orig.ident"
TK91 <- RenameIdents(TK91, "ATAC" = "E11")
TK91 <- AddMetaData(TK91, TK91@active.ident, "Age")
Idents(TK91) <- "orig.ident"
TK91 <- RenameIdents(TK91, "ATAC" = "TK91")
TK91 <- AddMetaData(TK91, TK91@active.ident, "Sample")
Idents(TK91) <- "Age"

load(file = "scATAC/Files/TK43.Robj")
Idents(TK43) <- "orig.ident"
TK43 <- RenameIdents(TK43, "ATAC" = "E12")
TK43 <- AddMetaData(TK43, TK43@active.ident, "Age")
Idents(TK43) <- "orig.ident"
TK43 <- RenameIdents(TK43, "ATAC" = "TK43")
TK43 <- AddMetaData(TK43, TK43@active.ident, "Sample")
Idents(TK43) <- "Age"

load(file = "scATAC/Files/TK47.Robj")
Idents(TK47) <- "orig.ident"
TK47 <- RenameIdents(TK47, "ATAC" = "E12")
TK47 <- AddMetaData(TK47, TK47@active.ident, "Age")
Idents(TK47) <- "orig.ident"
TK47 <- RenameIdents(TK47, "ATAC" = "TK47")
TK47 <- AddMetaData(TK47, TK47@active.ident, "Sample")
Idents(TK47) <- "Age"

load(file = "scATAC/Files/TK73.Robj")
Idents(TK73) <- "orig.ident"
TK73 <- RenameIdents(TK73, "ATAC" = "E13")
TK73 <- AddMetaData(TK73, TK73@active.ident, "Age")
Idents(TK73) <- "orig.ident"
TK73 <- RenameIdents(TK73, "ATAC" = "TK73")
TK73 <- AddMetaData(TK73, TK73@active.ident, "Sample")
Idents(TK73) <- "Age"

load(file = "scATAC/Files/TK95.Robj")
Idents(TK95) <- "orig.ident"
TK95 <- RenameIdents(TK95, "ATAC" = "E13")
TK95 <- AddMetaData(TK95, TK95@active.ident, "Age")
Idents(TK95) <- "orig.ident"
TK95 <- RenameIdents(TK95, "ATAC" = "TK95")
TK95 <- AddMetaData(TK95, TK95@active.ident, "Sample")
Idents(TK95) <- "Age"

load(file = "scATAC/Files/TK74.Robj")
Idents(TK74) <- "orig.ident"
TK74 <- RenameIdents(TK74, "ATAC" = "E14")
TK74 <- AddMetaData(TK74, TK74@active.ident, "Age")
Idents(TK74) <- "orig.ident"
TK74 <- RenameIdents(TK74, "ATAC" = "TK74")
TK74 <- AddMetaData(TK74, TK74@active.ident, "Sample")
Idents(TK74) <- "Age"

load(file = "scATAC/Files/TK75.Robj")
Idents(TK75) <- "orig.ident"
TK75 <- RenameIdents(TK75, "ATAC" = "E14")
TK75 <- AddMetaData(TK75, TK75@active.ident, "Age")
Idents(TK75) <- "orig.ident"
TK75 <- RenameIdents(TK75, "ATAC" = "TK75")
TK75 <- AddMetaData(TK75, TK75@active.ident, "Sample")
Idents(TK75) <- "Age"

load(file = "scATAC/Files/TKC4.Robj")
Idents(TKC4) <- "orig.ident"
TKC4 <- RenameIdents(TKC4, "ATAC" = "E14")
TKC4 <- AddMetaData(TKC4, TKC4@active.ident, "Age")
Idents(TKC4) <- "orig.ident"
TKC4 <- RenameIdents(TKC4, "ATAC" = "TKC4")
TKC4 <- AddMetaData(TKC4, TKC4@active.ident, "Sample")
Idents(TKC4) <- "Age"

####Fragment####
DefaultAssay(TK72) <- 'peaks'
DefaultAssay(TK91) <- 'peaks'
DefaultAssay(TK43) <- 'peaks'
DefaultAssay(TK47) <- 'peaks'
DefaultAssay(TK73) <- 'peaks'
DefaultAssay(TK95) <- 'peaks'
DefaultAssay(TK74) <- 'peaks'
DefaultAssay(TK75) <- 'peaks'
DefaultAssay(TKC4) <- 'peaks'

frags <- Fragments(TK72)  # get list of fragment objects
Fragments(TK72) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK72/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK72) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK91)  # get list of fragment objects
Fragments(TK91) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK91/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK91) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK43)  # get list of fragment objects
Fragments(TK43) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK43/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK43) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK47)  # get list of fragment objects
Fragments(TK47) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK47/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK47) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK73)  # get list of fragment objects
Fragments(TK73) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK73/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK73) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK95)  # get list of fragment objects
Fragments(TK95) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK95/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK95) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK74)  # get list of fragment objects
Fragments(TK74) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK74/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK74) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK75)  # get list of fragment objects
Fragments(TK75) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK75/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK75) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TKC4)  # get list of fragment objects
Fragments(TKC4) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TKC4/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TKC4) <- frags  # assign update list of fragment objects back to the assay

####Intersection####
# Find the union of peaks between the two datasets
# Extract peak regions from both datasets
DefaultAssay(TK72) <- 'peaks'
DefaultAssay(TK91) <- 'peaks'
DefaultAssay(TK43) <- 'peaks'
DefaultAssay(TK47) <- 'peaks'
DefaultAssay(TK73) <- 'peaks'
DefaultAssay(TK95) <- 'peaks'
DefaultAssay(TK74) <- 'peaks'
DefaultAssay(TK75) <- 'peaks'
DefaultAssay(TKC4) <- 'peaks'

peaks1 <- GetAssayData(TK72, assay = "peaks", slot = "ranges")
peaks2 <- GetAssayData(TK91, assay = "peaks", slot = "ranges")
peaks3 <- GetAssayData(TK43, assay = "peaks", slot = "ranges")
peaks4 <- GetAssayData(TK47, assay = "peaks", slot = "ranges")
peaks5 <- GetAssayData(TK73, assay = "peaks", slot = "ranges")
peaks6 <- GetAssayData(TK95, assay = "peaks", slot = "ranges")
peaks7 <- GetAssayData(TK74, assay = "peaks", slot = "ranges")
peaks8 <- GetAssayData(TK75, assay = "peaks", slot = "ranges")
peaks9 <- GetAssayData(TKC4, assay = "peaks", slot = "ranges")

# Find the intersection of peaks
combined.peaks <- reduce(x = c(peaks1, peaks2, peaks3, peaks4, peaks5, peaks6, peaks7, peaks8, peaks9))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# Subset both objects to include only common peaks
DefaultAssay(TK72) <- 'peaks'
DefaultAssay(TK91) <- 'peaks'
DefaultAssay(TK43) <- 'peaks'
DefaultAssay(TK47) <- 'peaks'
DefaultAssay(TK73) <- 'peaks'
DefaultAssay(TK95) <- 'peaks'
DefaultAssay(TK74) <- 'peaks'
DefaultAssay(TK75) <- 'peaks'
DefaultAssay(TKC4) <- 'peaks'

TK72.counts <- FeatureMatrix(
  fragments = Fragments(TK72),
  features = combined.peaks,
  cells = colnames(TK72))

TK91.counts <- FeatureMatrix(
  fragments = Fragments(TK91),
  features = combined.peaks,
  cells = colnames(TK91))

TK43.counts <- FeatureMatrix(
  fragments = Fragments(TK43),
  features = combined.peaks,
  cells = colnames(TK43))

TK47.counts <- FeatureMatrix(
  fragments = Fragments(TK47),
  features = combined.peaks,
  cells = colnames(TK47))

TK73.counts <- FeatureMatrix(
  fragments = Fragments(TK73),
  features = combined.peaks,
  cells = colnames(TK73))

TK95.counts <- FeatureMatrix(
  fragments = Fragments(TK95),
  features = combined.peaks,
  cells = colnames(TK95))

TK74.counts <- FeatureMatrix(
  fragments = Fragments(TK74),
  features = combined.peaks,
  cells = colnames(TK74))

TK75.counts <- FeatureMatrix(
  fragments = Fragments(TK75),
  features = combined.peaks,
  cells = colnames(TK75))

TKC4.counts <- FeatureMatrix(
  fragments = Fragments(TKC4),
  features = combined.peaks,
  cells = colnames(TKC4))

TK72[['peakunion']] <- CreateAssayObject(counts = TK72.counts)
TK91[['peakunion']] <- CreateAssayObject(counts = TK91.counts)
TK43[['peakunion']] <- CreateAssayObject(counts = TK43.counts)
TK47[['peakunion']] <- CreateAssayObject(counts = TK47.counts)
TK73[['peakunion']] <- CreateAssayObject(counts = TK73.counts)
TK95[['peakunion']] <- CreateAssayObject(counts = TK95.counts)
TK74[['peakunion']] <- CreateAssayObject(counts = TK74.counts)
TK75[['peakunion']] <- CreateAssayObject(counts = TK75.counts)
TKC4[['peakunion']] <- CreateAssayObject(counts = TKC4.counts)

DefaultAssay(TK72) <- 'peakunion'
DefaultAssay(TK91) <- 'peakunion'
DefaultAssay(TK43) <- 'peakunion'
DefaultAssay(TK47) <- 'peakunion'
DefaultAssay(TK73) <- 'peakunion'
DefaultAssay(TK95) <- 'peakunion'
DefaultAssay(TK74) <- 'peakunion'
DefaultAssay(TK75) <- 'peakunion'
DefaultAssay(TKC4) <- 'peakunion'

####Merging####
# Merge the datasets
pattern_atac <- merge(TK72, y = c(TK91, TK43, TK47, TK73, TK95, TK74, TK75, TKC4),
                       add.cell.ids = c("TK72", "TK91", "TK43", "TK47", "TK73", "TK95", "TK74", "TK75", "TKC4"))
rm(TK72)
rm(TK91)
rm(TK43)
rm(TK47)
rm(TK73)
rm(TK95)
rm(TK74)
rm(TK75)
rm(TKC4)

save(pattern_atac, file = "/media/thomaskim/Data/scATAC/Files/pattern_atac.Robj")

# Run TF-IDF normalization
pattern_atac <- RunTFIDF(pattern_atac)

# Find variable features
pattern_atac <- FindTopFeatures(pattern_atac, min.cutoff = "q0")

# Run dimensional reduction (LSI)
pattern_atac <- RunSVD(
  object = pattern_atac,
  assay = 'peakunion', #peaks
  reduction.key = 'LSI_',
  reduction.name = 'lsi')

pattern_atac <- RunHarmony(object = pattern_atac,
                            assay.use = 'peaks',
                            group.by.vars = 'Sample',
                            reduction.use = 'lsi', project.dim = F)

# Run UMAP
DepthCor(pattern_atac)
pattern_atac <- RunUMAP(pattern_atac, reduction = 'harmony', dims = 2:30)

# Find neighbors and clusters
pattern_atac <- FindNeighbors(pattern_atac, reduction = 'harmony', dims = 2:30)
#pattern_atac <- FindClusters(pattern_atac, verbose = FALSE)

# Plot UMAP
DimPlot(pattern_atac, label = , group.by = "Sample")
DimPlot(pattern_atac, label = , group.by = "Age")


save(pattern_atac, file = "/media/thomaskim/Data/scATAC/Files/pattern_atac.Robj")