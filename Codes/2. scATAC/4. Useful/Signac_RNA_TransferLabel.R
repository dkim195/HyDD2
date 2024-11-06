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
setwd("/media/thomaskim/Data/")

####Loading####
#ATAC
load(file = "scATAC/Files/TK43.Robj")
Idents(TK43) <- "orig.ident"
TK43 <- RenameIdents(TK43, "ATAC" = "E12")
TK43 <- AddMetaData(TK43, TK43@active.ident, "Age")
TK43 <- AddMetaData(TK43, TK43@active.ident, "Cluster_Pass2")
Idents(TK43) <- "Cluster_Pass2"

#RNA
load(file = "scRNA/Robj/E12_Atlas_Merge_V2.Robj")
E12_Atlas <- UpdateSeuratObject(E12_Atlas)
E12_Atlas@meta.data$old.ident <- NULL
E12_Atlas@meta.data$Cluster_Pass1 <- NULL
E12_Atlas@meta.data$Cluster <- NULL
E12_Atlas@meta.data$SCT_snn_res.0.4 <- NULL
E12_Atlas@meta.data$seurat_clusters <- NULL
Idents(E12_Atlas) <- "Age"
E12_Atlas <- RenameIdents(E12_Atlas, "E12" = "Atlas", "E12_Rep1" = "Atlas", "E12_Rep2" = "Atlas",
                          "E12_Atlas1" = "Atlas", "E12_Atlas2" = "Atlas", "E12_Atlas3" = "Atlas",
                          "E12_Atlas4" = "Atlas", "E12_Atlas5" = "Atlas", "E12_Atlas6" = "Atlas",
                          "E12_Atlas7" = "Atlas", "E12_Atlas8" = "Atlas", "E12_Atlas9" = "Atlas" )
E12_Atlas <- AddMetaData(E12_Atlas, E12_Atlas@active.ident, "Reference")
Idents(E12_Atlas) <- "Cluster_Pass2"

E12_Atlas <- RenameIdents(E12_Atlas, "Prethalamus (Sp8)" = "Prethalamus",
                          "Prethalamus (Sp9)" = "Prethalamus", "Prethalamus (Sst)" = "Prethalamus",
                          "PVH_SON (Cartpt)" = "PVH_SON", "PVH_SON (Otp)" = "PVH_SON")
E12_Atlas <- AddMetaData(E12_Atlas, E12_Atlas@active.ident, "Cluster_Pass2")

options(future.globals.maxSize = 8000 * 1024^2)
E12_Atlas <- E12_Atlas[, sample(colnames(E12_Atlas), size = 72930, replace=F)] #70000 seems to work
DefaultAssay(TK43) <- 'RNA'
DefaultAssay(E12_Atlas) <- 'RNA'

E12_Atlas <- FindVariableFeatures(object = E12_Atlas,  nfeatures = 1000)

transfer.anchors <- FindTransferAnchors(
  reference = E12_Atlas, query = TK43,
  normalization.method = "LogNormalize",
  reduction = "pcaproject", dims = 1:30)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = E12_Atlas$Cluster_Pass2,
  weight.reduction = TK43[['lsi']],
  dims = 2:30)

TK43 <- AddMetaData(object = TK43, metadata = predicted.labels)

DimPlot(object = TK43, label = TRUE, group.by = "predicted.id") + NoLegend()

FeaturePlot(
  object = TK43,
  features = c("Sp8","Nr4a2","Pitx2","Hmx2","Hmx3"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3)

