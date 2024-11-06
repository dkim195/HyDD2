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
library(stringr)
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 200000 * 1024^2) #200 gb ram

setwd("/faststorage/project/Hypothalamus")


####Mutant####
#Dlx1/2
load(file = "data/scRNA/Dlx1_2_Final.Robj")
#colnames(Mutant) <- sub("^TK[0-9]+_", "", colnames(Mutant))
Mutant@meta.data <- Mutant@meta.data[, c("orig.ident", "nCount_RNA","nFeature_RNA",
                                         "Cluster_Pass2", "Genotype"), drop = FALSE]

# Consolidate all 'Prethalamus_' types to 'Prethalamus'
Mutant@meta.data$Cluster_Pass2 <- str_replace(Mutant@meta.data$Cluster_Pass2, "Prethalamus.+", "Prethalamus")

# Consolidate all 'PVH_SON_' types to 'PVH_SON'
Mutant@meta.data$Cluster_Pass2 <- str_replace(Mutant@meta.data$Cluster_Pass2, "PVH_SON.+", "PVH_SON")

#For mutant, make another meta.data section combining clustes and genotype
Mutant$Clusters <- paste(Mutant$Genotype, Mutant$Cluster_Pass2, sep = "_")

#Cluster names should not have Remove parentheses, spaces, commas, underscores, and hyphens.
Mutant@meta.data$Clusters <- gsub("[() ,_-]", "", Mutant@meta.data$Clusters)

MuDataSeurat::WriteH5AD(Mutant, "SCENIC/RNA/OUTPUT/DLX/Mutant_Dlx.h5ad", assay="RNA")
write.csv(Embeddings(Mutant, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/DLX/Mutant_Dlx_cell_embeddings.csv")

rm(Mutant)
gc()

#Lhx2
load(file = "data/scRNA/Lhx2_Final.Robj")
#colnames(Mutant) <- sub("^TK[0-9]+_", "", colnames(Mutant))
Mutant@meta.data <- Mutant@meta.data[, c("orig.ident", "nCount_RNA","nFeature_RNA",
                                         "Cluster_Pass2", "Genotype"), drop = FALSE]

# Consolidate all 'Prethalamus_' types to 'Prethalamus'
Mutant@meta.data$Cluster_Pass2 <- str_replace(Mutant@meta.data$Cluster_Pass2, "Prethalamus.+", "Prethalamus")

# Consolidate all 'PVH_SON_' types to 'PVH_SON'
Mutant@meta.data$Cluster_Pass2 <- str_replace(Mutant@meta.data$Cluster_Pass2, "PVH_SON.+", "PVH_SON")

#For mutant, make another meta.data section combining clustes and genotype
Mutant$Clusters <- paste(Mutant$Genotype, Mutant$Cluster_Pass2, sep = "_")

#Cluster names should not have Remove parentheses, spaces, commas, underscores, and hyphens.
Mutant@meta.data$Clusters <- gsub("[() ,_-]", "", Mutant@meta.data$Clusters)

MuDataSeurat::WriteH5AD(Mutant, "SCENIC/RNA/OUTPUT/LHX2/Mutant_Lhx2.h5ad", assay="RNA")
write.csv(Embeddings(Mutant, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/LHX2/Mutant_Lhx2_cell_embeddings.csv")

rm(Mutant)
gc()

#Isl1
load(file = "data/scRNA/Isl1_Final.Robj")
#colnames(Mutant) <- sub("^TK[0-9]+_", "", colnames(Mutant))
Mutant@meta.data <- Mutant@meta.data[, c("orig.ident", "nCount_RNA","nFeature_RNA",
                                         "Cluster_Pass2", "Genotype"), drop = FALSE]

# Consolidate all 'Prethalamus_' types to 'Prethalamus'
Mutant@meta.data$Cluster_Pass2 <- str_replace(Mutant@meta.data$Cluster_Pass2, "Prethalamus.+", "Prethalamus")

# Consolidate all 'PVH_SON_' types to 'PVH_SON'
Mutant@meta.data$Cluster_Pass2 <- str_replace(Mutant@meta.data$Cluster_Pass2, "PVH_SON.+", "PVH_SON")

#For mutant, make another meta.data section combining clustes and genotype
Mutant$Clusters <- paste(Mutant$Genotype, Mutant$Cluster_Pass2, sep = "_")

#Cluster names should not have Remove parentheses, spaces, commas, underscores, and hyphens.
Mutant@meta.data$Clusters <- gsub("[() ,_-]", "", Mutant@meta.data$Clusters)

MuDataSeurat::WriteH5AD(Mutant, "SCENIC/RNA/OUTPUT/ISL1/Mutant_Isl1.h5ad", assay="RNA")
write.csv(Embeddings(Mutant, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/ISL1/Mutant_Isl1_cell_embeddings.csv")

rm(Mutant)
gc()

#Foxd1
load(file = "data/scRNA/Foxd1_Final.Robj")
#colnames(Mutant) <- sub("^TK[0-9]+_", "", colnames(Mutant))
Idents(Mutant) <- "Cluster_Pass2"
Mutant <- RenameIdents(Mutant, "Avp_Gal_Oxt (PVN_SON)" = "PVN_SON", "Bsx (PMN)"= "PMN", 
                    "Calb2 (SMN)" = "SMN", "Cck (MMN)" = "MMN",
                    "Check (POA_SCN)" = "POA_SCN", "Ghrh_Gal (PMN)"= "PMN", 
                    "Npvf_Grp_Hcrt (LH)" = "LH", "Npy_Pomc_Agrp (ARC)" = "ARC_VMH",
                    "NPC (Astro)" = "NPC (Glial)", 
                    "NPC (Epen)" = "NPC (Glial)", "NPC (Glial)" = "NPC (Glial)",
                    "NPC (Oligo)" = "NPC (Glial)", "NPC (Tanycyte)" = "NPC (Glial)",
                    "Nts_Tac1_Cck (MMN)" = "MMN", 
                    "Penk (POA_SCN)"= "POA_SCN", "Pmch_Trh (LH)"= "LH",
                    "Rora_Rorb (POA_SCN)"= "POA_SCN", "Tac2 (ARC_VMH)"= "ARC_VMH")
Mutant <- AddMetaData(Mutant, Mutant@active.ident, "Cluster_Pass2")

Mutant@meta.data <- Mutant@meta.data[, c("orig.ident", "nCount_RNA","nFeature_RNA",
                                         "Cluster_Pass2", "Genotype"), drop = FALSE]

#For mutant, make another meta.data section combining clustes and genotype
Mutant$Clusters <- paste(Mutant$Genotype, Mutant$Cluster_Pass2, sep = "_")

#Cluster names should not have Remove parentheses, spaces, commas, underscores, and hyphens.
Mutant@meta.data$Clusters <- gsub("[() ,_-]", "", Mutant@meta.data$Clusters)

MuDataSeurat::WriteH5AD(Mutant, "SCENIC/RNA/OUTPUT/FOXD1/Mutant_Foxd1.h5ad", assay="RNA")
write.csv(Embeddings(Mutant, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/FOXD1/Mutant_Foxd1_cell_embeddings.csv")

rm(Mutant)
gc()