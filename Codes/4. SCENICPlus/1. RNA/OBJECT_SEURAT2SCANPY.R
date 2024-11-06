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
library(MuDataSeurat)
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 200000 * 1024^2) #200 gb ram

setwd("/faststorage/project/Hypothalamus")

####Pattern####
load(file = "data/scRNA/pattern.Robj")
#colnames(pattern) <- sub("^TK[0-9]+_", "", colnames(pattern))
pattern@meta.data <- pattern@meta.data[, c("orig.ident", "nCount_RNA","nFeature_RNA",
                                         "Cluster_Pass2", "Genotype","Age_Sum"), drop = FALSE]
MuDataSeurat::WriteH5AD(pattern, "SCENIC/RNA/OUTPUT/pattern.h5ad", assay="RNA")
write.csv(Embeddings(pattern, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/pattern_cell_embeddings.csv")
rm(pattern)
gc()



####Neurogenesis####
load(file = "data/scRNA/neurogenesis.Robj")
#colnames(neurogenesis) <- sub("^TK[0-9]+_", "", colnames(neurogenesis))
neurogenesis@meta.data <- neurogenesis@meta.data[, c("orig.ident", "nCount_RNA","nFeature_RNA",
                                         "Cluster_Pass2", "Genotype","Age_Sum"), drop = FALSE]
MuDataSeurat::WriteH5AD(neurogenesis, "SCENIC/RNA/OUTPUT/neurogenesis.h5ad", assay="RNA")
write.csv(Embeddings(neurogenesis, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/neurogenesis_cell_embeddings.csv")
rm(neurogenesis)
gc()




####Mutant####
#Dlx1/2
load(file = "data/scRNA/Dlx1_2_Final.Robj")
#colnames(Mutant) <- sub("^TK[0-9]+_", "", colnames(Mutant))
Mutant@meta.data <- Mutant@meta.data[, c("orig.ident", "nCount_RNA","nFeature_RNA",
                                         "Cluster_Pass2", "Genotype"), drop = FALSE]
Idents(Mutant) <- "Genotype"
CTRL <- subset(Mutant, idents = "Control")
CKO <- subset(Mutant, idents = "Mutant")

MuDataSeurat::WriteH5AD(CTRL, "SCENIC/RNA/OUTPUT/DLX/CTRL_DLX.h5ad", assay="RNA")
write.csv(Embeddings(CTRL, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/DLX/CTRL_DLX_cell_embeddings.csv")

MuDataSeurat::WriteH5AD(CKO, "SCENIC/RNA/OUTPUT/DLX/CKO_DLX.h5ad", assay="RNA")
write.csv(Embeddings(CKO, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/DLX/CKO_DLX_cell_embeddings.csv")
rm(Mutant)
rm(CTRL)
rm(CKO)
gc()

#Lhx2
load(file = "data/scRNA/Lhx2_Final.Robj")
#colnames(Mutant) <- sub("^TK[0-9]+_", "", colnames(Mutant))
Mutant@meta.data <- Mutant@meta.data[, c("orig.ident", "nCount_RNA","nFeature_RNA",
                                         "Cluster_Pass2", "Genotype"), drop = FALSE]
Idents(Mutant) <- "Genotype"
CTRL <- subset(Mutant, idents = "Control")
CKO <- subset(Mutant, idents = "Mutant")

MuDataSeurat::WriteH5AD(CTRL, "SCENIC/RNA/OUTPUT/Lhx2/CTRL_Lhx2.h5ad", assay="RNA")
write.csv(Embeddings(CTRL, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/Lhx2/CTRL_Lhx2_cell_embeddings.csv")

MuDataSeurat::WriteH5AD(CKO, "SCENIC/RNA/OUTPUT/Lhx2/CKO_Lhx2.h5ad", assay="RNA")
write.csv(Embeddings(CKO, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/Lhx2/CKO_Lhx2_cell_embeddings.csv")
rm(Mutant)
rm(CTRL)
rm(CKO)
gc()

#Isl1
load(file = "data/scRNA/Isl1_Final.Robj")
#colnames(Mutant) <- sub("^TK[0-9]+_", "", colnames(Mutant))
Mutant@meta.data <- Mutant@meta.data[, c("orig.ident", "nCount_RNA","nFeature_RNA",
                                         "Cluster_Pass2", "Genotype"), drop = FALSE]
Idents(Mutant) <- "Genotype"
CTRL <- subset(Mutant, idents = "Control")
CKO <- subset(Mutant, idents = "Mutant")

MuDataSeurat::WriteH5AD(CTRL, "SCENIC/RNA/OUTPUT/Isl1/CTRL_Isl1.h5ad", assay="RNA")
write.csv(Embeddings(CTRL, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/Isl1/CTRL_Isl1_cell_embeddings.csv")

MuDataSeurat::WriteH5AD(CKO, "SCENIC/RNA/OUTPUT/Isl1/CKO_Isl1.h5ad", assay="RNA")
write.csv(Embeddings(CKO, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/Isl1/CKO_Isl1_cell_embeddings.csv")
rm(Mutant)
rm(CTRL)
rm(CKO)
gc()

#Foxd1
load(file = "data/scRNA/Foxd1_Final.Robj")
#colnames(Mutant) <- sub("^TK[0-9]+_", "", colnames(Mutant))
Mutant@meta.data <- Mutant@meta.data[, c("orig.ident", "nCount_RNA","nFeature_RNA",
                                         "Cluster_Pass2", "Genotype"), drop = FALSE]
Idents(Mutant) <- "Genotype"
CTRL <- subset(Mutant, idents = "Control")
CKO <- subset(Mutant, idents = "Mutant")

MuDataSeurat::WriteH5AD(CTRL, "SCENIC/RNA/OUTPUT/Foxd1/CTRL_Foxd1.h5ad", assay="RNA")
write.csv(Embeddings(CTRL, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/Foxd1/CTRL_Foxd1_cell_embeddings.csv")

MuDataSeurat::WriteH5AD(CKO, "SCENIC/RNA/OUTPUT/Foxd1/CKO_Foxd1.h5ad", assay="RNA")
write.csv(Embeddings(CKO, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/Foxd1/CKO_Foxd1_cell_embeddings.csv")
rm(Mutant)
rm(CTRL)
rm(CKO)
gc()