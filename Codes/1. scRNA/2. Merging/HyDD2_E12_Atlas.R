library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)

####HyDD1_Prep####
#E12 atlas
setwd("D:/")

#TK15
TK15 <- Read10X(data.dir = "scRNA/TK15-NkxHET1/raw_gene_bc_matrices/mm10/")
colnames(TK15) = paste0("TK15_", colnames(TK15))
TK15 <- CreateSeuratObject(counts = TK15, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK15 <- subset(TK15, subset = nCount_RNA > 1000)  
TK15 <- RenameIdents(TK15, "TK15"= "E12_Atlas1")
TK15 <- AddMetaData(TK15, TK15@active.ident, "Age")
Idents(TK15) <- "Age"
TK15[["percent.mt"]] <- PercentageFeatureSet(TK15, pattern = "^mt-")
TK15[["percent.RPS"]] <- PercentageFeatureSet(TK15, pattern = "^Rps")
TK15[["percent.RPL"]] <- PercentageFeatureSet(TK15, pattern = "^Rpl")
TK15 <- subset(TK15, subset = percent.mt <50)
TK15 <- subset(TK15, subset = percent.RPS <25)

#TK16
TK16 <- Read10X(data.dir = "scRNA/TK16-NkxHET2/raw_gene_bc_matrices/mm10/")
colnames(TK16) = paste0("TK16_", colnames(TK16))
TK16 <- CreateSeuratObject(counts = TK16, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK16 <- subset(TK16, subset = nCount_RNA > 1000)  
TK16 <- RenameIdents(TK16, "TK16"= "E12_Atlas2")
TK16 <- AddMetaData(TK16, TK16@active.ident, "Age")
Idents(TK16) <- "Age"
TK16[["percent.mt"]] <- PercentageFeatureSet(TK16, pattern = "^mt-")
TK16[["percent.RPS"]] <- PercentageFeatureSet(TK16, pattern = "^Rps")
TK16[["percent.RPL"]] <- PercentageFeatureSet(TK16, pattern = "^Rpl")
TK16 <- subset(TK16, subset = percent.mt <50)
TK16 <- subset(TK16, subset = percent.RPS <25)

#TK17
TK17 <- Read10X(data.dir = "scRNA/TK17-Ex3HET1/raw_gene_bc_matrices/mm10/")
colnames(TK17) = paste0("TK17_", colnames(TK17))
TK17 <- CreateSeuratObject(counts = TK17, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK17 <- subset(TK17, subset = nCount_RNA > 1000)  
TK17 <- RenameIdents(TK17, "TK17"= "E12_Atlas3")
TK17 <- AddMetaData(TK17, TK17@active.ident, "Age")
Idents(TK17) <- "Age"
TK17[["percent.mt"]] <- PercentageFeatureSet(TK17, pattern = "^mt-")
TK17[["percent.RPS"]] <- PercentageFeatureSet(TK17, pattern = "^Rps")
TK17[["percent.RPL"]] <- PercentageFeatureSet(TK17, pattern = "^Rpl")
TK17 <- subset(TK17, subset = percent.mt <50)
TK17 <- subset(TK17, subset = percent.RPS <25)

#TK19
TK19 <- Read10X(data.dir = "scRNA/TK19/raw_gene_bc_matrices/mm10/")
colnames(TK19) = paste0("TK19_", colnames(TK19))
TK19 <- CreateSeuratObject(counts = TK19, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK19 <- subset(TK19, subset = nCount_RNA > 1000)  
TK19 <- RenameIdents(TK19, "TK19"= "E12_Atlas4")
TK19 <- AddMetaData(TK19, TK19@active.ident, "Age")
Idents(TK19) <- "Age"
TK19[["percent.mt"]] <- PercentageFeatureSet(TK19, pattern = "^mt-")
TK19[["percent.RPS"]] <- PercentageFeatureSet(TK19, pattern = "^Rps")
TK19[["percent.RPL"]] <- PercentageFeatureSet(TK19, pattern = "^Rpl")
TK19 <- subset(TK19, subset = percent.mt <50)
TK19 <- subset(TK19, subset = percent.RPS <25)

#TK22
TK22 <- Read10X(data.dir = "scRNA/TK22/raw_gene_bc_matrices/mm10/")
colnames(TK22) = paste0("TK22_", colnames(TK22))
TK22 <- CreateSeuratObject(counts = TK22, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK22 <- subset(TK22, subset = nCount_RNA > 1000)  
TK22 <- RenameIdents(TK22, "TK22"= "E12_Atlas5")
TK22 <- AddMetaData(TK22, TK22@active.ident, "Age")
Idents(TK22) <- "Age"
TK22[["percent.mt"]] <- PercentageFeatureSet(TK22, pattern = "^mt-")
TK22[["percent.RPS"]] <- PercentageFeatureSet(TK22, pattern = "^Rps")
TK22[["percent.RPL"]] <- PercentageFeatureSet(TK22, pattern = "^Rpl")
TK22 <- subset(TK22, subset = percent.mt <50)
TK22 <- subset(TK22, subset = percent.RPS <25)

#TK25
TK25 <- Read10X(data.dir = "scRNA/TK25/raw_gene_bc_matrices/mm10/")
colnames(TK25) = paste0("TK25_", colnames(TK25))
TK25 <- CreateSeuratObject(counts = TK25, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK25 <- subset(TK25, subset = nCount_RNA > 1000)  
TK25 <- RenameIdents(TK25, "TK25"= "E12_Atlas6")
TK25 <- AddMetaData(TK25, TK25@active.ident, "Age")
Idents(TK25) <- "Age"
TK25[["percent.mt"]] <- PercentageFeatureSet(TK25, pattern = "^mt-")
TK25[["percent.RPS"]] <- PercentageFeatureSet(TK25, pattern = "^Rps")
TK25[["percent.RPL"]] <- PercentageFeatureSet(TK25, pattern = "^Rpl")
TK25 <- subset(TK25, subset = percent.mt <50)
TK25 <- subset(TK25, subset = percent.RPS <25)

#TK32
TK32 <- Read10X(data.dir = "scRNA/TK32/raw_gene_bc_matrices/mm10/")
colnames(TK32) = paste0("TK32_", colnames(TK32))
TK32 <- CreateSeuratObject(counts = TK32, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK32 <- subset(TK32, subset = nCount_RNA > 1000)  
TK32 <- RenameIdents(TK32, "TK32"= "E12_Atlas7")
TK32 <- AddMetaData(TK32, TK32@active.ident, "Age")
Idents(TK32) <- "Age"
TK32[["percent.mt"]] <- PercentageFeatureSet(TK32, pattern = "^mt-")
TK32[["percent.RPS"]] <- PercentageFeatureSet(TK32, pattern = "^Rps")
TK32[["percent.RPL"]] <- PercentageFeatureSet(TK32, pattern = "^Rpl")
TK32 <- subset(TK32, subset = percent.mt <50)
TK32 <- subset(TK32, subset = percent.RPS <25)

#TK35
TK35 <- Read10X(data.dir = "scRNA/TK35/outs/raw_feature_bc_matrix/")
colnames(TK35) = paste0("TK35_", colnames(TK35))
TK35 <- CreateSeuratObject(counts = TK35, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK35 <- subset(TK35, subset = nCount_RNA > 1000)  
TK35 <- RenameIdents(TK35, "TK35"= "E12_Atlas8")
TK35 <- AddMetaData(TK35, TK35@active.ident, "Age")
Idents(TK35) <- "Age"
TK35[["percent.mt"]] <- PercentageFeatureSet(TK35, pattern = "^mt-")
TK35[["percent.RPS"]] <- PercentageFeatureSet(TK35, pattern = "^Rps")
TK35[["percent.RPL"]] <- PercentageFeatureSet(TK35, pattern = "^Rpl")
TK35 <- subset(TK35, subset = percent.mt <50)
TK35 <- subset(TK35, subset = percent.RPS <25)

#TK36
TK36 <- Read10X(data.dir = "scRNA/TK36/outs/raw_feature_bc_matrix/")
colnames(TK36) = paste0("TK36_", colnames(TK36))
TK36 <- CreateSeuratObject(counts = TK36, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK36 <- subset(TK36, subset = nCount_RNA > 1000)  
TK36 <- RenameIdents(TK36, "TK36"= "E12_Atlas9")
TK36 <- AddMetaData(TK36, TK36@active.ident, "Age")
Idents(TK36) <- "Age"
TK36[["percent.mt"]] <- PercentageFeatureSet(TK36, pattern = "^mt-")
TK36[["percent.RPS"]] <- PercentageFeatureSet(TK36, pattern = "^Rps")
TK36[["percent.RPL"]] <- PercentageFeatureSet(TK36, pattern = "^Rpl")
TK36 <- subset(TK36, subset = percent.mt <50)
TK36 <- subset(TK36, subset = percent.RPS <25)

#TK15 1203 cells 1950 UMI 900 genes
#TK16 195 cells 1600 UMI 800 genes
#TK17 48 cells 2300 UMI 1000 genes
#TK19 14335 cells 4400 UMI 1800 genes
#TK22 8762 cells 4200 UMI 1800 genes
#TK25 7111 cells, 4620 UMI, 1800 genes
#TK32 8823 cells 3900 UMI 1600 genes
#TK35 11622 cells 1900 UMI 900 genes
#TK36 14484 cells 2100 UMI 950 genes



#merge all and clean up before merging to HyDD2 E12
E12_Atlas <- merge(x = TK15, y = list(TK16, TK17, TK19, TK22, TK25, TK32, TK35, TK36))
E12_Atlas[["percent.sex"]] <- PercentageFeatureSet(E12_Atlas,
                                                   features = c("Xist","Malat1","Tsix"))
VlnPlot(E12_Atlas, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                                "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)                 
save(E12_Atlas, file = "Robj/E12_Atlas_Rest.Robj")

#process
E12_Atlas <- SCTransform(E12_Atlas, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
E12_Atlas<- RunPCA(E12_Atlas, npcs = 50, ndims.print = NA, verbose = F)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
m.s.genes$mouseGene
m.g2m.genes$mouseGene
E12_Atlas <- CellCycleScoring(E12_Atlas, s.features = m.s.genes$mouseGene,
                              g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
E12_Atlas<- RunHarmony(E12_Atlas, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E12_Atlas<- RunUMAP(E12_Atlas, reduction = "harmony", dims = 1:20,
                    n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(E12_Atlas, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(E12_Atlas, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
DimPlot(E12_Atlas, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(E12_Atlas, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
        split.by = "Phase", ncol = 2) + NoLegend()

E12_Atlas <- FindNeighbors(E12_Atlas, dims = 1:20, reduction = "harmony")
E12_Atlas <- FindClusters(E12_Atlas, resolution = 2) #SCT_snn_res.0.8
DimPlot(E12_Atlas, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
markers <- FindAllMarkers(E12_Atlas, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/E12_Atlas-Cluster-DEG-initial.csv")
save(E12_Atlas, file = "Robj/E12_Atlas_Rest.Robj")


#0 PVH&SON
#1 PRETHAL
#2 NPC
#3 G2M NPC/NPC
#4 PRETHAL
#5 JUNK
#6 NEURLA PRO ASCL1
#7 NEURLA PRO NEUROG2
#8 PRETHAL
#9 SMN/MMN progenitor
#10 G2M NPC/NPC
#11 G2M NPC/NPC
#12 SMN
#13 JUNK
#14 NEURAL PRO?
#15 JUNK
#16 TUBERAL
#17 Thalamic progenitor
#18 NPC
#19 SMN
#20 MMN
#21 JUNK
#22 PVH&SON
#23 JUNK
#24 TUBERAL PMN
#25 NPC
#26 PRETHAL
#27 Thalamic eminence
#28 Immune cells
#29 JUNK
#30 G2M NPC/NPC
#31 JUNK
#32 PRETHAL
#33 PRETHAL
#34 SMN/MMN progenitor
#35 SMN
#36 NEURLA PRO ASCL1
#37 JUNK
#38 PVH&SON
#39 PVH&SON
#40 JUNK
#41 JUNK
#42 NEURLA PRO NEUROG2
#43 JUNK-IMMUNE
#44 ZLI? 



#first clean up
E12_Atlas <- RenameIdents(E12_Atlas, "0" = "PVH/SON", "1" = "Prethalamus", "2" = "NPC", "3" = "G2M NPC/NPC",
                          "4" = "Prethalamus", "5" = "Junk",  "6" = "Neural Pro (Ascl1)",
                          "7" = "Neural Pro (Neurog2)", "8" =  "Prethalamus",
                          "9" = "Neural Pro (SMN/MMN)", "10" =  "G2M NPC/NPC",
                          "11" = "G2M NPC/NPC", "12" = "SMN", "13" = "Junk", "14" = "Neural Pro", "15" = "Junk",
                          "16" = "Tuberal (ARC/VMH)", "17" = "Thalamic progenitor",
                          "18" = "NPC", "19" = "SMN",
                          "20" ="MMN", "21" = "Junk", "22" = "PVH/SON",
                          "23" = "Junk", "24" = "Tuberal (PMN)", "25" ="NPC", "26" = "Prethalamus",
                          "27" = "Thalamic eminence", "28" = "Immune cells", "29" = "Junk", "30" = "G2M NPC/NPC",
                          "31" = "Junk", "32" = "Prethalamus", "33" = "Prethalamus",
                          "34" = "Neural Pro (SMN/MMN)",
                          "35" = "SMN", "36" = "Neural Pro (Ascl1)",
                          "37" = "Junk", "38" = "PVH/SON",
                          "39" = "PVH/SON", "40" = "Junk",
                          "41" = "Junk", "42" = "Neural Pro (Neurog2)", 
                          "43" = "Immune cells", "44" = "ZLI" )
Group <- sort(levels(E12_Atlas@active.ident))
E12_Atlas <- AddMetaData(E12_Atlas, E12_Atlas@active.ident, "Cluster_Pass1")
E12_Atlas@active.ident <- factor(x = E12_Atlas@active.ident, levels = Group)
E12_Atlas@meta.data$Cluster_Pass1 <- factor(x = E12_Atlas@meta.data$Cluster_Pass1, levels = Group)
DimPlot(E12_Atlas, reduction = "umap", label = F, pt.size = 0.5)
save(E12_Atlas, file = "Robj/E12_Atlas_Rest.Robj")

#save Immune cell for later
E12_Atlas <- subset(E12_Atlas, idents = c("PVH/SON", "Prethalamus", "NPC","G2M NPC/NPC",
                                          "Neural Pro (Ascl1)","Neural Pro (Neurog2)",
                                          "Neural Pro (SMN/MMN)","SMN",  "Neural Pro",
                                          "Tuberal (ARC/VMH)", "Thalamic progenitor",
                                          "MMN", "Tuberal (PMN)",  "Thalamic eminence"))
E12_Atlas <- SCTransform(E12_Atlas, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
E12_Atlas<- RunPCA(E12_Atlas, npcs = 50, ndims.print = NA, verbose = F)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
E12_Atlas <- CellCycleScoring(E12_Atlas, s.features = m.s.genes$mouseGene,
                              g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
E12_Atlas<- RunHarmony(E12_Atlas, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E12_Atlas<- RunUMAP(E12_Atlas, reduction = "harmony", dims = 1:20,
                    n.neighbours = 10L, min.dist = 0.01, spread = 1)
Idents(E12_Atlas) <- "Cluster_Pass1"
DimPlot(E12_Atlas, reduction = "umap", label = T, )
DimPlot(E12_Atlas, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(E12_Atlas, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
save(E12_Atlas, file = "Robj/E12_Atlas_Rest_V2.Robj")


####E12 HyDD2####
#E12 merge file 
load(file = "Robj/E12_Merge.Robj")
DimPlot(E12, reduction = "umap", label = F, pt.size = 0.5)
E12 <- subset(E12, idents = c("SMN", "Prethalamus", "G2M NPC/NPC", 
                    "PVH/SON", "Check (Junk?)", "MMN",  "NPC",
                    "Neural Pro (Ascl1)", "Neural Pro (Neurog2)",
                    "Tuberal (PMN)", "Prethalamus (ID)", "Foxg1/PVH/SON",
                     "Thalamic eminence", "Tuberal (ARC/VMH)",
                    "Tuberal (PMN-LH)", 
                    "Neural Pro (SMN/MMN)",  "Foxg1 (NPC/Thalamic eminence?)"))
DimPlot(E12, reduction = "umap", label = F, pt.size = 0.5)
E12@meta.data$Age <- E12@meta.data$orig.ident

####Merging####
E12_Atlas <- merge(x = E12, y = list(E12_Atlas))
E12_Atlas@meta.data$old.ident <- NULL
E12_Atlas@meta.data$SCT_snn_res.2 <- NULL
E12_Atlas@meta.data$seurat_clusters <- NULL
VlnPlot(E12_Atlas, features = c("nCount_RNA","nFeature_RNA"), pt.size = 0, group.by = "Age", ncol = 1)

save(E12_Atlas, file = "Robj/E12_Atlas_Merge.Robj")
load(file = "Robj/E12_Atlas_Merge.Robj")

E12_Atlas <- SCTransform(E12_Atlas,
                         vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
E12_Atlas<- RunPCA(E12_Atlas, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(E12_Atlas, ndims = 50, reduction = "pca" )
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
E12_Atlas <- CellCycleScoring(E12_Atlas, s.features = m.s.genes$mouseGene,
                              g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
E12_Atlas<- RunHarmony(E12_Atlas, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E12_Atlas<- RunUMAP(E12_Atlas, reduction = "harmony", dims = 1:28,
                    n.neighbours = 10L, min.dist = 0.01, spread = 1)

Idents(E12_Atlas) <- "Cluster_Pass1"
DimPlot(E12_Atlas, reduction = "umap", label = T) + NoLegend() + NoAxes()

####Check_Cluster####
E12_Atlas <- FindNeighbors(E12_Atlas, dims = 1:28, reduction = "harmony")
E12_Atlas <- FindClusters(E12_Atlas, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(E12_Atlas, reduction = "umap", label = T, pt.size = 0.5) + NoLegend() + NoAxes()

DimPlot(E12_Atlas, reduction = "umap", label = T, group.by = "Cluster_Pass1") + NoLegend() + NoAxes()

markers <- FindAllMarkers(E12_Atlas, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/E12_Atlas-Cluster-DEG-initial.csv")
save(E12_Atlas, file = "Robj/E12_Atlas_Merge.Robj")


####Final####
#Some Foxg1 but not great deal,
#might capture some Bnst
#0 NPC
#1 SMN
#2 Neural Progenitor
#3 Neural Progenitor (Neurog2)
#4 Prethalamus (Sp8)
#5 Prethalamus (Sp9)
#6 PVH & SON (Otp)
#7 Junk
#8 MMN
#9 PMN
#10 Neural Progenitor (Ascl1)
#11 PVH&SON (Cartpt)
#12 AntID_ID_TT
#13 Tuberal (ARC/VMH)
#14 Thalamic Progenitor
#15 Thalamic Eminence
#16 Prethalamus (Sst)
#17 Junk?
#18 Junk? Ascl1 G2M NPC might be doublet
#19 LH (Pmch)


E12_Atlas <- RenameIdents(E12_Atlas, "0" = "NPC", "1" = "SMN", "2" = "Neural Pro", "3" = "Neural Pro (Neurog2)",
                          "4" = "Prethalamus (Sp8)", "5" = "Prethalamus (Sp9)", "6" = "PVH_SON (Otp)",
                          "7" = "Junk", "8" = "MMN", "9" = "PMN", "10" = "Neural Pro (Ascl1)", 
                          "11" = "PVH_SON (Cartpt)", "12" = "AntID_ID_TT", "13" = "Tuberal (ARC_VMH)",
                          "14" = "Thalamic progenitor", "15" = "Thalamic eminence", "16" = "Prethalamus (Sst)",
                          "17" = "Junk", "18" = "Junk", "19" = "LH (Pmch)")
Group <- sort(levels(E12_Atlas@active.ident))
E12_Atlas <- AddMetaData(E12_Atlas, E12_Atlas@active.ident, "Cluster_Pass2")
E12_Atlas@active.ident <- factor(x = E12_Atlas@active.ident, levels = Group)
E12_Atlas@meta.data$Cluster_Pass1 <- factor(x = E12_Atlas@meta.data$Cluster_Pass1, levels = Group)
DimPlot(E12_Atlas, reduction = "umap", label = F, pt.size = 0.5)


#remove junk
E12_Atlas <- subset(E12_Atlas, idents = c("NPC", "SMN", "Neural Pro", "Neural Pro (Neurog2)",
                                          "Prethalamus (Sp8)", "Prethalamus (Sp9)", "PVH_SON (Otp)",
                                          "MMN", "PMN", "Neural Pro (Ascl1)", 
                                          "PVH_SON (Cartpt)","AntID_ID_TT", "Tuberal (ARC_VMH)",
                                          "Thalamic progenitor", "Thalamic eminence","Prethalamus (Sst)",
                                          "LH (Pmch)"))
E12_Atlas <- SCTransform(E12_Atlas,
                         vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
E12_Atlas<- RunPCA(E12_Atlas, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(E12_Atlas, ndims = 50, reduction = "pca" )
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
E12_Atlas <- CellCycleScoring(E12_Atlas, s.features = m.s.genes$mouseGene,
                              g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
E12_Atlas<- RunHarmony(E12_Atlas, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E12_Atlas<- RunUMAP(E12_Atlas, reduction = "harmony", dims = 1:25,
                    n.neighbours = 10L, min.dist = 0.01, spread = 1)

Idents(E12_Atlas) <- "Cluster_Pass2"
DimPlot(E12_Atlas, reduction = "umap", label = T) + NoLegend() + NoAxes()

DimPlot(E12_Atlas, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(E12_Atlas, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
DimPlot(E12_Atlas, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(E12_Atlas, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
        split.by = "Phase", ncol = 2) + NoLegend()

Group <- sort(levels(E12_Atlas@active.ident))
E12_Atlas <- AddMetaData(E12_Atlas, E12_Atlas@active.ident, "Cluster_Pass2")
E12_Atlas@active.ident <- factor(x = E12_Atlas@active.ident, levels = Group)
E12_Atlas@meta.data$Cluster_Pass2 <- factor(x = E12_Atlas@meta.data$Cluster_Pass2, levels = Group)
DimPlot(E12_Atlas, reduction = "umap", label = T, pt.size = 0.5)

save(E12_Atlas, file = "Robj/E12_Atlas_Merge.Robj")

markers <- FindAllMarkers(E12_Atlas, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E12_Atlas-Cluster-DEG-initial.csv")

#test the weird junk island
E12_Atlas <- FindNeighbors(E12_Atlas, dims = 1:25, reduction = "harmony")
E12_Atlas <- FindClusters(E12_Atlas, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(E12_Atlas, reduction = "umap", label = T, pt.size = 0.5) + NoLegend() + NoAxes()

markers <- FindAllMarkers(E12_Atlas, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E12_Atlas-Cluster-DEG-initial_Clusterpass2.csv")


#Remove cluster 11
#Remove Thalamic progenitor c14
#Remove Thalamic eminence c13
E12_Atlas <- subset(E12_Atlas, idents = c("0","1","2","3","4","5","6","7","8",
                                          "9","10","12","15","16","17","18"))
E12_Atlas <- SCTransform(E12_Atlas,
                         vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
E12_Atlas<- RunPCA(E12_Atlas, npcs = 50, ndims.print = NA, verbose = F)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
E12_Atlas <- CellCycleScoring(E12_Atlas, s.features = m.s.genes$mouseGene,
                              g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
E12_Atlas<- RunHarmony(E12_Atlas, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E12_Atlas<- RunUMAP(E12_Atlas, reduction = "harmony", dims = 1:22,
                    n.neighbours = 10L, min.dist = 0.01, spread = 1)
Idents(E12_Atlas) <- "Cluster_Pass2"
DimPlot(E12_Atlas, reduction = "umap", label = T) + NoLegend() + NoAxes()


#FINAL CLEAN UP
#Remove Thalamic eminence and progenitor
table(E12_Atlas@active.ident)
E12_Atlas <- subset(E12_Atlas, idents = c("NPC", "SMN", "Neural Pro", "Neural Pro (Neurog2)",
                                          "Prethalamus (Sp8)", "Prethalamus (Sp9)", "PVH_SON (Otp)",
                                          "MMN", "PMN", "Neural Pro (Ascl1)", 
                                          "PVH_SON (Cartpt)","AntID_ID_TT", "Tuberal (ARC_VMH)",
                                          "Prethalamus (Sst)","LH (Pmch)"))
Group <- sort(levels(E12_Atlas@active.ident))
E12_Atlas <- AddMetaData(E12_Atlas, E12_Atlas@active.ident, "Cluster_Pass2")
E12_Atlas@active.ident <- factor(x = E12_Atlas@active.ident, levels = Group)
E12_Atlas@meta.data$Cluster_Pass2 <- factor(x = E12_Atlas@meta.data$Cluster_Pass2, levels = Group)
DimPlot(E12_Atlas, reduction = "umap", label = T, pt.size = 0.5)

#remove extra thalamic eminence and progenitor
Idents(E12_Atlas) <- "Cluster_Pass1"
E12_Atlas@meta.data$Test <- paste(Idents(E12_Atlas), E12_Atlas$Phase, sep = "_")
table(E12_Atlas@meta.data$Test)
Idents(E12_Atlas) <- "Test"
E12_Atlas <- subset(E12_Atlas, idents = c("MMN_G1", "MMN_G2M", "MMN_S", "NA_G1", 
                                          "NA_G2M", "NA_S", "Neural Pro (Ascl1)_G1",  "Neural Pro (Ascl1)_G2M", 
                                          "Neural Pro (Ascl1)_S", "Neural Pro (Neurog2)_G1",
                                          "Neural Pro (Neurog2)_G2M", "Neural Pro (Neurog2)_S", 
                                          "Neural Pro_G1", "Neural Pro_G2M", "Neural Pro_S", "NPC_G1", 
                                          "NPC_G2M", "NPC_S", "SMN_G1", "SMN_G2M", "SMN_S"))
Idents(E12_Atlas) <- "Cluster_Pass2"
Group <- sort(levels(E12_Atlas@active.ident))
E12_Atlas <- AddMetaData(E12_Atlas, E12_Atlas@active.ident, "Cluster_Pass2")
E12_Atlas@active.ident <- factor(x = E12_Atlas@active.ident, levels = Group)
E12_Atlas@meta.data$Cluster_Pass2 <- factor(x = E12_Atlas@meta.data$Cluster_Pass2, levels = Group)
DimPlot(E12_Atlas, reduction = "umap", label = T, pt.size = 0.5)
E12_Atlas@meta.data$Test <- NULL

FeaturePlot(E12_Atlas, c("Tbr1","Lhx9","Lhx2","Rspo3","Neurod6","Neurog2"))

save(E12_Atlas, file = "Robj/E12_Atlas_Merge_V2.Robj")

markers <- FindAllMarkers(E12_Atlas, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E12_Atlas-Cluster-DEG.csv")




Cells1 <-WhichCells(E12_Atlas, idents = "Thalamic progenitor")
DimPlot(E12_Atlas, reduction = "umap", label = F, 
        pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
  ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf")) 










####Plot####
genes <- c("Cdkn1c","Gadd45g","Pkib","Dlk1","Hes6","Ass1") #AH Neural Progenitor
genes <- c("D930028M14Rik", "Onecut2", "Sncg", "Snhg11", "Ina", "Cartpt") #Ant ID
genes <- c("Pomc", "Chchd10", "S100a10", "Six6","Rbp1","Sox14") #ARC
genes <- c("Sst","Isl1","Six6","Otp","Dlk1","Cited1") #DMH
genes <- c("Lhx8", "Snhg11","Cadm1","Pcsk1n","Mapt","Gabra2","Lhx6") #ID/TT
genes <- c("Apoe","Lgals1","Igfbp7","Sepp1","Serpinh1","Tyrobp") #Immune cell
genes <- c("Pmch","Lhx9","Chd3os") #LH?
genes <- c("Dlx2","Lhx8","Dlx6os1","Sp9","Arx","Foxg1") #MGE
genes <- c("Lhx1os","Pcp4","Foxb1","Lhx1","Nhlh2","Lhx5") #MMN
genes <- c("Ube2c","Ccnb1","Hmgb2","Spc25","Cks2","Neurog2") #G2M NPC
genes <- c("Pitx2","Gm28050","Nrn1","Uncx","Tmem163","Nhlh2","Pax7") #Pax7 TT
genes <- c("Hmx2","Bub3","Cited1","Isl1","Hmx3") #PMN
genes <- c("Six6","Lhx8","Gad2","Ina","Dlx5","Syt1") #POA
genes <- c("2810417H13Rik", "Fabp7", "Rrm2","Hes5", "Hmgb2","H2afz") #Prethalamic NPC
genes <- c("Meis2","Isl1","Pax6","Sp8", "Gad2", "Arx") #Prethalamus
genes <- c("Otp", "Cartpt","Tbca","Arxes2","Gm2694","Fezf1") #PVH & SON
genes <- c("Tcf7l2","Lhx1os","Gata3","Tal2","Gata2","Lhx1") #Rim domain
genes <- c("Zic1","Dlx5","Arx","Gad2","Sp9","Calb1") #SCN
genes <- c("Pitx2","Barhl1","Nr4a2","Lmx1a","Chgb","Foxa1") #SMN
genes <- c("Neurog2","Gadd45g","Cdkn1c","Nhlh2","Ppp1r14a","Neurog1") #SMN/MMN progenitor
genes <- c("Tbr1","Rspo3","Calb2","Lhx9","Nrn1","Neurod6") #Thalamic eminence
genes <- c("Tcf7l2","Lhx2","Shox2","Lhx9","Neurog2","Gbx1") #Thalamic progenitor
genes <- c("Chchd10","Gm2694","Fezf1","Snca","Sox14","Plp1") #VMH
genes <- c("Dbx1","Id3","Wnt8b","Tpbg","Hes1","Fgfbp3") #ZLI
genes <- c("Fabp7","Hes5", "Hes1", "Fabp7","Id3","Hmgb2") #Hypo NPC

p <- FeaturePlot(E12_Atlas, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = F, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)



####Plan####

#Check cluster

#Use this as a reference data base on seurat integration

#option 1
#subset based 1. subset(object, subset = nFeature_RNA > 2000) and 2. subset(object, subset = nCount_RNA > 2000)
#rerun SCTransform, PCA, Harmony, UMAP
#downsample to match mutant as  subset(object, cells = sample(Cells(object), 1000)

#option 2
#just downsample to match mutant as  subset(object, cells = sample(Cells(object), 1000)
