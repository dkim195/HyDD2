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
path <- getwd()
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 512000 * 1024^2) #512 gb ram
setwd("/faststorage/project/Hypothalamus/data")

####Saving####
args = commandArgs(trailingOnly=TRUE)
load("scATAC/OBJECT.Robj")

ct <- OBJECT@assays$peaks@counts
row.names(ct)<-sub("-", ":", row.names(ct))
#Signac chr14-89896606-89897308 -> cistopic chr14:89896606-89897308
ct <- as.data.frame(ct)
ct$rownames <- rownames(ct)
write_feather(ct, sink=paste0(path, '/OBJECT_count_matrix.feather'))
write.csv(Cells(OBJECT), file = "OBJECT_cellID_obs.csv")
write.csv(Embeddings(OBJECT, reduction = "umap"), file = "OBJECT_cell_embeddings.csv")
write.table(OBJECT@meta.data, file = "OBJECT_cell_annotation.tsv",
            sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)

