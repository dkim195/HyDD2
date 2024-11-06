
library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)

####E16####
setwd("D:/")
load(file = "Robj/E16_Merge.Robj")

clusters <- levels(E16$Cluster_Pass1)
fix(clusters)
E16<- subset(E16, idents = c("ARC/VMH", "Avp/Gal/Oxt (PVN/SON?)", "Bsx/Cck (PMN?)", "Cck (MMN)", 
                                          "Check (Ant ID?)", "Check (PMN?)", "Check (SMN/PVN/SON?)", "Check (SMN?)", 
                                          "Gal (PMN)", "Glial Pro?", "Glial Pro? (Tanycyte?)", "Hcrt/Npvf", 
                                           "MMN", "Nts/Tac1/Cck (MMN)", "Penk (POA/SCN?)", 
                                          "PMN", "Pomc/Npy/Sst (Arc)", "PreThal/ID", "PreThal/ID (Glial?)", 
                                          "Rora/Rorb (POA/SCN?)", "Tac1 (SMN)", "Trh/Pnoc (POA/SCN?)"))
E16<- SCTransform(E16, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
E16<- RunPCA(E16, npcs = 50, ndims.print = NA, verbose = F)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
E16<- CellCycleScoring(E16, s.features = m.s.genes$mouseGene,
                              g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
E16<- RunHarmony(E16, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E16<- RunUMAP(E16, reduction = "harmony", dims = 1:20,
                    n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(E16, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#Second clustering

E16 <- FindNeighbors(E16, dims = 1:20, reduction = "harmony")
E16 <- FindClusters(E16, resolution = 2) #SCT_snn_res.0.8
DimPlot(E16, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(E16, file = "Robj/E16_Atlas.Robj")

markers <- FindAllMarkers(E16, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/DEG-E16-Atlas.csv")

cluster <- levels(E16@active.ident)
fix(cluster)
E16 <- RenameIdents(E16, "0" = "ARC_VMH", "1" = "PMN", "2" = "Ghrh_Gal (PMN)",
                    "3" = "Cck (MMN)", "4" = "PreThal", "5" = "Sst_Npy", "6" = "Npy_Pomc_Agrp (ARC)",
                    "7" = "Calb2 (SMN)", "8" = "MMN", "9" = "Bsx (PMN)",
                    "10" = "Penk (POA_SCN)", "11" = "Neural Pro (Ascl1)", 
                    "12" = "NPC (Glial)", "13" = "SMN", "14" = "NPC (Epen)",
                    "15"="MMN", "16"="Nts_Tac1_Cck (MMN)", "17"="MMN", "18"="NPC (Tanycyte)",
                    "19"="AntID_ID", "20"="NPC (Epen)", "21" = "Calb2 (SMN)", "22"="SMN", 
                    "23"="Pmch_Trh (LH)", "24"="ARC_VMH", "25"="Rora_Rorb (POA_SCN)",
                    "26"="Check (POA_SCN)", "27"="PreThal_ID", "28"="Neural Pro (Neurog2)",
                    "29"="NPC (Oligo)", "30"="NPC (Astro)", "31"="NPC (Glial)",
                    "32"="NPC (Oligo)", "33"="NPC (Glial)", 
                    "34"="Avp_Gal_Oxt (PVN_SON)", "35"="NPC (Epen)", "36"="Npvf_Grp_Hcrt (LH)",
                    "37"="MMN", "38"="NPC (Astro)", "39"="Tac2 (ARC_VMH)", "40"="NPC (Oligo)")
Group <- sort(levels(E16@active.ident))
E16 <- AddMetaData(E16, E16@active.ident, "Cluster_Pass1")
E16@active.ident <- factor(x = E16@active.ident, levels = Group)
E16@meta.data$Cluster_Pass1 <- factor(x = E16@meta.data$Cluster_Pass1, levels = Group)
DimPlot(E16, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

save(E16, file = "Robj/E16_Atlas.Robj")
markers <- FindAllMarkers(E16, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/DEG-E16-Atlas.csv")


