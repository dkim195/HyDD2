library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("D:/")
load(file = "Robj/neurogenesis.Robj")
Idents(neurogenesis) <- "Age_Sum"
neurogenesis <- subset(neurogenesis, idents = "P8")

Idents(neurogenesis) <- "Cluster_Pass2"
neurogenesis <- subset(neurogenesis, idents = c("ARC (Npy_Agrp)", "ARC (Pomc)", "ARC_VMH", "ARC_VMH (Tac2_PMN)", 
                                                "DMH-LH (Grp_Cck)", "DMH-LH (Npvf)", "DMH (Npw_PMN)", 
                                                "Isl1 cells", "LH (Hcrt_Oxt)", "LH (Pmch_Trh)", "MMN", "MMN (Cck)", 
                                                "MMN (Npy_Cck)", "MMN (Nts_Tac1_Cck)", "PMN (Ghrh_Gal_Cited1)", 
                                                "PMN (Ghrh_Gal_Oxt)", "POA_SCN", "PreThal",
                                                "PreThal_ID", "PVN_SON (Avp_Gal_Oxt)", 
                                                "SCN (Rorb)", "SMN", "SMN (Calb2)", "SMN (Tac2)", "Sst_Npy", 
                                                "TMN_PH (Hdc)"))
setwd("E:/")
load(file = "Robj/TKC5_TKC6.Robj")
TKC5_TKC6 <- subset(TKC5_TKC6, idents = "Neuron")

P8 <- merge(x = neurogenesis, y = list(TKC5_TKC6))
rm(neurogenesis)
rm(TKC5_TKC6)
head(P8@meta.data)
P8@meta.data$percent.mt <- NULL
P8@meta.data$percent.RPS <- NULL
P8@meta.data$percent.RPL <- NULL
P8@meta.data$percent.sex <- NULL
P8@meta.data$S.Score <- NULL
P8@meta.data$G2M.Score <- NULL
P8@meta.data$Phase <- NULL
P8@meta.data$old.ident <- NULL
P8@meta.data$SCT_snn_res.0.4 <- NULL
P8@meta.data$SCT_snn_res.2 <- NULL
P8@meta.data$seurat_clusters <- NULL

P8[["percent.sex"]] <- PercentageFeatureSet(P8,features = c("Xist","Malat1","Tsix"))
P8[["percent.mt"]] <- PercentageFeatureSet(P8, pattern = "^mt-")
P8[["percent.RPS"]] <- PercentageFeatureSet(P8, pattern = "^Rps")
P8[["percent.RPL"]] <- PercentageFeatureSet(P8, pattern = "^Rpl")

P8<- SCTransform(P8, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
P8<- RunPCA(P8, npcs = 50, ndims.print = NA, verbose = F)

P8<- RunHarmony(P8, group.by.vars = c("orig.ident"), assay.use = "SCT", plot.convergence = T)
P8<- RunUMAP(P8, reduction = "harmony", dims = 1:20,
                  n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(P8, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

P8 <- FindNeighbors(P8, dims = 1:20, reduction = "harmony")
P8 <- FindClusters(P8, resolution = 2) #SCT_snn_res.0.8
DimPlot(P8, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
markers <- FindAllMarkers(P8, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/P8_Neuron-Cluster-DEG-initial.csv")
save(P8, file = "Robj/P8_Neuron.Robj")

#0 Junk
#1 Junk
#2 PreThal
#3 ARC (Pomc)
#4 VMH
#5 SCN (Vip)
#6 SMN (Calb2_Tac1)
#7 MMN (Calb1_Cck)
#8 Junk
#9 PreThal
#10 Junk
#11 DMH-LH (Grp_Cck)
#12 Junk
#13 SCN (Rorb)
#14 MMN (Cartpt_Cck)
#15 MMN (Pcp4_Cck)
#16 MMN (Nts_Cartpt)
#17 MMN (Cck)
#18 PMN (Ghrh_Gal_Oxt)
#19 Junk
#20 Sst (TMN-PH?)
#21 ARC (Npy_Agrp)
#22 MMN (Nts_Tac1_Cck)
#23 PVN_SON (Avp_Gal_Oxt)
#24 Junk
#25 VMH (Tac1)
#26 reThal
#27 LH (Pmch_Trh)
#28 Junk
#29 Junk
#30 LH (Pmch_Tac2)
#31 MMN (Pcp4_Cck)
#32 Tac2
#33 SMN (Nts)
#34 Sst (ARC-TMN?)
#35 PMN (Ghrh_Gal_Oxt)
#36 DMH (Npw_PMN)
#37 Junk
#38 LH (Hcrt)
#39 TMN_PH (Hdc)
#40 Junk
#41 SMN (Calb2)
#42 Sst_Npy
#43 SCN (Avp)
#44 Junk
#45 PVN_SON (Avp_Gal_Oxt)
#46 Junk
#47 ARC (Npw_Agrp)


names <- levels(P8@active.ident)
P8 <- RenameIdents(P8, c("0" = "Junk", "1" = "Junk", "2" = "PreThal", "3" = "ARC (Pomc)", "4" = "VMH",
                         "5" = "SCN (Vip)", "6" = "SMN (Calb2_Tac1)", "7" = "MMN (Calb1_Cck)", "8" = "Junk",
                         "9" = "PreThal", "10" = "Junk", "11" = "DMH-LH (Grp_Cck)", 
                         "12" = "Junk", "13" = "SCN (Rorb)",
                         "14" = "MMN (Cartpt_Cck)", "15" = "MMN (Pcp4_Cck)",
                         "16" = "MMN (Nts_Cartpt)", "17" = "MMN (Cck)", "18" = "PMN (Ghrh_Gal_Oxt)",
                         "19" = "Junk", "20" = "Sst (TMN-PH?)",
                         "21" = "ARC (Npy_Agrp)", "22" = "MMN (Nts_Tac1_Cck)", 
                         "23" = "PVN_SON (Avp_Gal_Oxt)", "24" = "Junk", "25" = "VMH (Tac1)", 
                         "26" = "PreThal", "27" = "LH (Pmch_Trh)", "28" = "Junk",
                         "29" = "Junk", "30" = "LH (Pmch_Tac2)", "31" = "MMN (Pcp4_Cck)",
                         "32" = "Tac2", "33" = "SMN (Nts)", 
                         "34" = "Sst (ARC-TMN?)", "35" = "PMN (Ghrh_Gal_Oxt)", "36" = "DMH (Npw_PMN)",
                         "37" = "Junk", "38" = "LH (Hcrt)", "39" = "TMN_PH (Hdc)",
                         "40" = "Junk", "41" = "SMN (Calb2)", "42" = "Sst_Npy",
                         "43" = "SCN (Avp)", "44" = "Junk", 
                         "45" = "PVN_SON (Avp_Gal_Oxt)", "46" = "Junk", "47" = "ARC (Npw_Agrp)"))
Group <- sort(levels(P8@active.ident))
P8 <- AddMetaData(P8, P8@active.ident, "Cluster_Pass3")
P8@active.ident <- factor(x = P8@active.ident, levels = Group)
P8@meta.data$Cluster_Pass3<- factor(x = P8@meta.data$Cluster_Pass3, levels = Group)
DimPlot(P8, reduction = "umap", label = F, pt.size = 0.5)

P8 <- subset(P8, idents = c("ARC (Npw_Agrp)", "ARC (Npy_Agrp)", "ARC (Pomc)", "DMH-LH (Grp_Cck)", 
                            "DMH (Npw_PMN)",  "LH (Hcrt)", "LH (Pmch_Tac2)", "LH (Pmch_Trh)", 
                            "MMN (Calb1_Cck)", "MMN (Cartpt_Cck)", "MMN (Cck)", "MMN (Nts_Cartpt)", 
                            "MMN (Nts_Tac1_Cck)", "MMN (Pcp4_Cck)", "PMN (Ghrh_Gal_Oxt)", 
                            "PreThal", "PVN_SON (Avp_Gal_Oxt)", "SCN (Avp)", "SCN (Rorb)", 
                            "SCN (Vip)", "SMN (Calb2)", "SMN (Calb2_Tac1)", "SMN (Nts)", 
                            "Sst (ARC-TMN?)", "Sst (TMN-PH?)", "Sst_Npy", "Tac2", "TMN_PH (Hdc)", 
                            "VMH", "VMH (Tac1)"))
P8<- SCTransform(P8, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
P8<- RunPCA(P8, npcs = 50, ndims.print = NA, verbose = F)
P8<- RunHarmony(P8, group.by.vars = c("orig.ident"), assay.use = "SCT", plot.convergence = T)
P8<- RunUMAP(P8, reduction = "harmony", dims = 1:20,
             n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(P8, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(P8, file = "Robj/P8_Neuron.Robj")

