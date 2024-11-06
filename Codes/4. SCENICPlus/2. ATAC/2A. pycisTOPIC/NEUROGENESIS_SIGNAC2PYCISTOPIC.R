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

#neurogenesis
neurogenesis_path <- file.path(path, "neurogenesis")
load("scATAC/neurogenesis_atac_annotation_v3b.Robj")

#reduce NPC (Glia) as there are already way too many cells in that clusters
Idents(neurogenesis_atac) <- "Cluster_Pass2"
# Step 1: Identify cells labeled as "NPC (Glia)"
npc_glia_cells <- which(neurogenesis_atac@meta.data$Cluster_Pass2 == "NPC (Glia)")
# Step 2: Randomly sample 10% of these cells
set.seed(42) # for reproducibility
npc_glia_sample <- sample(npc_glia_cells, length(npc_glia_cells) * 0.1)
# Step 3: Include the sampled "NPC (Glia)" cells and all other cells
all_other_cells <- WhichCells(neurogenesis_atac, expression = Cluster_Pass2 != "NPC (Glia)")
cells_to_keep <- c(npc_glia_sample, all_other_cells)
# Create a new Seurat object with the selected cells
neurogenesis_atac <- subset(neurogenesis_atac, cells = cells_to_keep)

#rename
neurogenesis_atac <- RenameIdents(neurogenesis_atac, 
                                  "PreThal_AntID" = "PreThal",
                                  "PreThal_AntID_ZI" = "PreThal")
neurogenesis_atac <- AddMetaData(neurogenesis_atac, neurogenesis_atac@active.ident, "Cluster_Pass2")
#Cluster names should not have Remove parentheses, spaces, commas, underscores, and hyphens.
neurogenesis_atac@meta.data$Cluster_Pass2 <- gsub("[() ,_-]", "", neurogenesis_atac@meta.data$Cluster_Pass2)

#reduce number
neurogenesis_atac <- neurogenesis_atac[, sample(colnames(neurogenesis_atac), size = 80000, replace=F)] #size = cell no similar to mutant

gc()

ct <- neurogenesis_atac@assays$peaks@counts
row.names(ct)<-sub("-", ":", row.names(ct))
#Signac chr14-89896606-89897308 -> cistopic chr14:89896606-89897308
ct <- as.data.frame(ct)
ct$rownames <- rownames(ct)

write_feather(ct, sink=file.path(neurogenesis_path, "neurogenesis_count_matrix.feather"))
write.csv(Cells(neurogenesis_atac), file = file.path(neurogenesis_path, "neurogenesis_cellID_obs.csv"))
write.csv(Embeddings(neurogenesis_atac, reduction = "umap"), file = file.path(neurogenesis_path, "neurogenesis_cell_embeddings.csv"))
write.table(neurogenesis_atac@meta.data, file = file.path(neurogenesis_path, "neurogenesis_cell_annotation.tsv"),
            sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
rm(neurogenesis_atac)
rm(ct)
gc()
