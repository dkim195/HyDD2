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
plan("multicore", workers = 20)
options(future.globals.maxSize = 512000 * 1024^2) #512 gb ram
setwd("/faststorage/project/Hypothalamus/data")
path <- "/faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT"

####Saving####
args = commandArgs(trailingOnly=TRUE)

#PATTERN
PATTERN_path <- file.path(path, "PATTERN")
load("scATAC/pattern_atac_annotation_v3b.Robj")
#pattern_atac <- pattern_atac[, sample(colnames(pattern_atac), size = 50000, replace=F)] #size = cell no similar to mutant

ct <- pattern_atac@assays$peaks@counts
row.names(ct)<-sub("-", ":", row.names(ct))
#Signac chr14-89896606-89897308 -> cistopic chr14:89896606-89897308
ct <- as.data.frame(ct)
ct$rownames <- rownames(ct)

Idents(pattern_atac) <- "Cluster_Pass2"

#Cluster names should not have Remove parentheses, spaces, commas, underscores, and hyphens.
pattern_atac@meta.data$Cluster_Pass2 <- gsub("[() ,_-]", "", pattern_atac@meta.data$Cluster_Pass2)

write_feather(ct, sink=file.path(PATTERN_path, "PATTERN_count_matrix.feather"))
write.csv(Cells(pattern_atac), file = file.path(PATTERN_path, "PATTERN_cellID_obs.csv"))
write.csv(Embeddings(pattern_atac, reduction = "umap"), file = file.path(PATTERN_path, "PATTERN_cell_embeddings.csv"))
write.table(pattern_atac@meta.data, file = file.path(PATTERN_path, "PATTERN_cell_annotation.tsv"),
            sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
rm(pattern_atac)
rm(ct)
gc()
