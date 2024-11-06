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


####neurogenesis####
#Dlx1/2
load("data/scRNA/neurogenesis.Robj")
#colnames(neurogenesis) <- sub("^TK[0-9]+_", "", colnames(neurogenesis))
Idents(neurogenesis) <- "Cluster_Pass2"
neurogenesis <- RenameIdents(neurogenesis,
                             "ARC (Npy_Agrp)" = "ARC_VMH",
                             "ARC (Pomc)" = "ARC_VMH_PMN",
                             "ARC_VMH (Tac2_PMN)" = "ARC_VMH_KISS1", 
                             "DMH-LH (Grp_Cck)" = "PMN", 
                             "DMH-LH (Npvf)" = "Npvf_Lhx9",
                             "Ependymal" = "Tanycyte_Ependymal",
                             "LH (Hcrt_Oxt)" = "Npvf_Lhx9",
                             "LH (Pmch_Trh)" = "ARC_VMH_PMN_PMCH",
                             "MMN (Cck)" = "MMN",
                             "MMN (Npy_Cck)" = "MMN",
                             "MMN (Nts_Tac1_Cck)" = "MMN",
                             "Neural Pro (Ascl1/Neurog2)" = "Neural Prog (Neurog2)",
                             "Neural Pro (Neurog2)" = "Neural Prog (Neurog2)",
                             "Oligodendrocytes (newly)" = "Oligodendrocyte", 
                             "PMN (Ghrh_Gal_Cited1)" = "PMN", 
                             "PMN (Ghrh_Gal_Oxt)" = "PMN", 
                             "POA_SCN" = "SCN?", 
                             "PreThal_ID" = "PreThal",
                             "PVN_SON (Avp_Gal_Oxt)" = "PVH_SON",
                             "SCN (Rorb)" = "SCN?",  
                             "SMN (Calb2)" = "SMN",
                             "SMN (Tac2)" = "SMN",
                             "Tanycyte" = "Tanycyte_Ependymal", 
                             "TMN_PH (Hdc)" = "MMN",
                             "NPC (Epen)" = "NPC (Glia)",
                             "NPC (Epen/Tan)" = "NPC (Glia)", 
                             "NPC (Glial)" = "NPC (Glia)",
                             "NPC (Glial_Astro)" = "NPC (Glia)", 
                             "NPC (Glial_Check)" = "NPC (Glia)",
                             "NPC (Glial_Oligo)" = "NPC (Glia)",
                             "Astrocytes" = "Astrocyte")
neurogenesis <- AddMetaData(neurogenesis, neurogenesis@active.ident, "Cluster_Pass2")

neurogenesis@meta.data <- neurogenesis@meta.data[, c("orig.ident", "nCount_RNA","nFeature_RNA",
                                         "Cluster_Pass2", "Age_Sum"), drop = FALSE]

#Cluster names should not have Remove parentheses, spaces, commas, underscores, and hyphens.
neurogenesis@meta.data$Cluster_Pass2 <- gsub("[() ,_-]", "", neurogenesis@meta.data$Cluster_Pass2)

MuDataSeurat::WriteH5AD(neurogenesis, "SCENIC/RNA/OUTPUT/neurogenesis/neurogenesis.h5ad", assay="RNA")
write.csv(Embeddings(neurogenesis, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/neurogenesis/neurogenesis_cell_embeddings.csv")

rm(neurogenesis)
gc()

