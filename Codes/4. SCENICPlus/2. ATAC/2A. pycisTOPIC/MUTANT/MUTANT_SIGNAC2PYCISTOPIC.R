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
library(cisTopic)
library(arrow)
library(stringr)
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 512000 * 1024^2) #512 gb ram
setwd("/faststorage/project/Hypothalamus/data")
path <- "/faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT"

####Saving####
args = commandArgs(trailingOnly=TRUE)

#DLX
dlx_path <- file.path(path, "DLX")
load("scATAC/Dlx_Mutant_adjusted_Final.Robj")
ct <- Mutant_adjusted@assays$peaks@counts
row.names(ct)<-sub("-", ":", row.names(ct))
#Signac chr14-89896606-89897308 -> cistopic chr14:89896606-89897308
ct <- as.data.frame(ct)
ct$rownames <- rownames(ct)

#For Mutant, make another meta.data section combining clustes and genotype
Mutant_adjusted$Clusters <- paste(Mutant_adjusted$Genotype, Mutant_adjusted$Cluster_Pass2, sep = "_")

#Cluster names should not have Remove parentheses, spaces, commas, underscores, and hyphens.
Mutant_adjusted@meta.data$Clusters <- gsub("[() ,_-]", "", Mutant_adjusted@meta.data$Clusters)

write_feather(ct, sink=file.path(dlx_path, "Mutant_Dlx_count_matrix.feather"))
write.csv(Cells(Mutant_adjusted), file = file.path(dlx_path, "Mutant_Dlx_cellID_obs.csv"))
write.csv(Embeddings(Mutant_adjusted, reduction = "umap"), file = file.path(dlx_path, "Mutant_Dlx_cell_embeddings.csv"))
write.table(Mutant_adjusted@meta.data, file = file.path(dlx_path, "Mutant_Dlx_cell_annotation.tsv"),
            sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
rm(Mutant_adjusted)
rm(ct)
gc()

#Lhx2
lhx2_path <- file.path(path, "LHX2")
load("scATAC/Lhx2_Mutant_adjusted_Final.Robj")
ct <- Mutant_adjusted@assays$peaks@counts
row.names(ct)<-sub("-", ":", row.names(ct))
#Signac chr14-89896606-89897308 -> cistopic chr14:89896606-89897308
ct <- as.data.frame(ct)
ct$rownames <- rownames(ct)

#For Mutant, make another meta.data section combining clustes and genotype
Mutant_adjusted$Clusters <- paste(Mutant_adjusted$Genotype, Mutant_adjusted$Cluster_Pass2, sep = "_")

#Cluster names should not have Remove parentheses, spaces, commas, underscores, and hyphens.
Mutant_adjusted@meta.data$Clusters <- gsub("[() ,_-]", "", Mutant_adjusted@meta.data$Clusters)

write_feather(ct, sink=file.path(lhx2_path, "Mutant_Lhx2_count_matrix.feather"))
write.csv(Cells(Mutant_adjusted), file = file.path(lhx2_path, "Mutant_Lhx2_cellID_obs.csv"))
write.csv(Embeddings(Mutant_adjusted, reduction = "umap"), file = file.path(lhx2_path, "Mutant_Lhx2_cell_embeddings.csv"))
write.table(Mutant_adjusted@meta.data, file = file.path(lhx2_path, "Mutant_Lhx2_cell_annotation.tsv"),
            sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
rm(Mutant_adjusted)
rm(ct)
gc()

#Isl1
isl1_path <- file.path(path, "ISL1")
load("scATAC/Isl1_Mutant_adjusted_Final.Robj")
ct <- Mutant_adjusted@assays$peaks@counts
row.names(ct)<-sub("-", ":", row.names(ct))
#Signac chr14-89896606-89897308 -> cistopic chr14:89896606-89897308
ct <- as.data.frame(ct)
ct$rownames <- rownames(ct)

#For Mutant, make another meta.data section combining clustes and genotype
Mutant_adjusted$Cluster_Pass2 <- Mutant_adjusted$predicted.id
Mutant_adjusted$Clusters <- paste(Mutant_adjusted$Genotype, Mutant_adjusted$Cluster_Pass2, sep = "_")

#Cluster names should not have Remove parentheses, spaces, commas, underscores, and hyphens.
Mutant_adjusted@meta.data$Clusters <- gsub("[() ,_-]", "", Mutant_adjusted@meta.data$Clusters)

write_feather(ct, sink=file.path(isl1_path, "Mutant_Isl1_count_matrix.feather"))
write.csv(Cells(Mutant_adjusted), file = file.path(isl1_path, "Mutant_Isl1_cellID_obs.csv"))
write.csv(Embeddings(Mutant_adjusted, reduction = "umap"), file = file.path(isl1_path, "Mutant_Isl1_cell_embeddings.csv"))
write.table(Mutant_adjusted@meta.data, file = file.path(isl1_path, "Mutant_Isl1_cell_annotation.tsv"),
            sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
rm(Mutant_adjusted)
rm(ct)
gc()

#Foxd1
foxd1_path <- file.path(path, "FOXD1")
load("scATAC/Foxd1_Mutant_adjusted_Final.Robj")
ct <- Mutant_adjusted@assays$peaks@counts
row.names(ct)<-sub("-", ":", row.names(ct))
#Signac chr14-89896606-89897308 -> cistopic chr14:89896606-89897308
ct <- as.data.frame(ct)
ct$rownames <- rownames(ct)

#For Mutant, make another meta.data section combining clustes and genotype
Mutant_adjusted$Cluster_Pass2 <- Mutant_adjusted$predicted.id
Mutant_adjusted$Clusters <- paste(Mutant_adjusted$Genotype, Mutant_adjusted$Cluster_Pass2, sep = "_")

#Cluster names should not have Remove parentheses, spaces, commas, underscores, and hyphens.
Mutant_adjusted@meta.data$Clusters <- gsub("[() ,_-]", "", Mutant_adjusted@meta.data$Clusters)

write_feather(ct, sink=file.path(foxd1_path, "Mutant_Foxd1_count_matrix.feather"))
write.csv(Cells(Mutant_adjusted), file = file.path(foxd1_path, "Mutant_Foxd1_cellID_obs.csv"))
write.csv(Embeddings(Mutant_adjusted, reduction = "umap"), file = file.path(foxd1_path, "Mutant_Foxd1_cell_embeddings.csv"))
write.table(Mutant_adjusted@meta.data, file = file.path(foxd1_path, "Mutant_Foxd1_cell_annotation.tsv"),
            sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
rm(Mutant_adjusted)
rm(ct)
gc()

