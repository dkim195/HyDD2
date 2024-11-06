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


####pattern####
#Dlx1/2
load("data/scRNA/pattern.Robj")
#colnames(pattern) <- sub("^TK[0-9]+_", "", colnames(pattern))
Idents(pattern) <- "Cluster_Pass2"
pattern <- RenameIdents(pattern, "Ant_ID" = "AntID_ID", 
                        "MMN (Prog)" = "Neural Pro (Neurog2) 2", 
                        "PVH_SON_AntVen" = "PVN_SON",
                        "PVH_SON" = "PVN_SON",
                        "SST (DMH?)" = "SST (DMH)",
                        "NPC (gliogenic?)" = "NPC (Gliogenic)",
                        "NSC"  =  "NSC_NPC", 
                        "PMN (Pro)" = "Tub_PMN (Prog)",
                        "Tub (Prog)"= "Tub_PMN (Prog)",
                        "Neural Pro (Neurog2)" = "Neural Pro (Neurog2) 1",
                        "PMN? (Gsx1)" = "PMN",
                        "PreThal (Pro)" = "PreThal", 
                        "PreThal (Sp8)" = "PreThal",
                        "PreThal (Sp9)" = "PreThal",
                        "PreThal_AntID" = "PreThal",
                        "Tub (ARC_VMH)" = "Tub (ARC_VMH)", 
                        "Tub (LH-Pmch)" = "Tub (ARC_VMH)", 
                        "Tub (LH)"= "Tub (ARC_VMH)", 
                        "Neural Pro (Ascl1)" = "Hypo_PreThal NPC",
                        "NPC (Ascl1)" = "Hypo_PreThal NPC",
                        "NPC (Neurog2)" = "Hypo_PreThal NPC")
pattern <- AddMetaData(pattern, pattern@active.ident, "Cluster_Pass2")

pattern@meta.data <- pattern@meta.data[, c("orig.ident", "nCount_RNA","nFeature_RNA",
                                         "Cluster_Pass2", "Age_Sum"), drop = FALSE]

#Cluster names should not have Remove parentheses, spaces, commas, underscores, and hyphens.
pattern@meta.data$Cluster_Pass2 <- gsub("[() ,_-]", "", pattern@meta.data$Cluster_Pass2)

MuDataSeurat::WriteH5AD(pattern, "SCENIC/RNA/OUTPUT/PATTERN/pattern.h5ad", assay="RNA")
write.csv(Embeddings(pattern, reduction = "umap"), file = "SCENIC/RNA/OUTPUT/PATTERN/pattern_cell_embeddings.csv")

rm(pattern)
gc()

