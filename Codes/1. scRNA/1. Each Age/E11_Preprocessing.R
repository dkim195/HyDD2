library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)

####loading####
#Processing E11
#fresh
setwd("D:/")
#Load HyDDV1, subset E11
load(file = "HyDDv1/HyDD_pattern.Robj")
Idents(Hypo2) <- "orig.ident"
Hypo2 <- subset(Hypo2, idents = "E11")
Hypo2@meta.data$S.Score <- NULL
Hypo2@meta.data$G2M.Score <- NULL
Hypo2@meta.data$Phase <- NULL
Hypo2@meta.data$SCT_snn_res.1.2 <- NULL
Hypo2@meta.data$seurat_clusters <- NULL
Hypo2@meta.data$SCT_snn_res.2 <- NULL
Hypo2@meta.data$old.ident <- NULL
Hypo2@meta.data$NEWCLUSTER <- NULL
Hypo2@meta.data$protocol <- NULL
Hypo2@meta.data$nGene <- NULL
Hypo2@meta.data$nUMI <- NULL
Hypo2 <- AddMetaData(Hypo2, Hypo2@active.ident, "Age")
Hypo2 <- subset(Hypo2, subset = nFeature_RNA > 500)
Hypo2 <- subset(Hypo2, subset = nCount_RNA > 1000)
QUALITY <- Hypo2@meta.data
write.csv(QUALITY, file = "CSV/E11_HyDD1_QC.csv")
save(Hypo2, file = "Robj/E11_HyDDV1.Robj")

#New v3.1 10x HyDD2
TK64<- Read10X_h5("scRNA/TK64/outs/raw_feature_bc_matrix.h5")
#adding column names to avoid barcode clash
colnames(TK64) = paste0("TK64_", colnames(TK64))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK64 <- CreateSeuratObject(counts = TK64, project = "E11",
                           min.cells = 5, min.features = 1000) 
TK64 <- subset(TK64, subset = nCount_RNA > 2000) 

###merge####
TK64 <- RenameIdents(TK64, "TK64"= "E11_Rep")
TK64 <- AddMetaData(TK64, TK64@active.ident, "Age")
head(TK64@meta.data)
Idents(TK64) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK64[["percent.mt"]] <- PercentageFeatureSet(TK64, pattern = "^mt-")
TK64[["percent.RPS"]] <- PercentageFeatureSet(TK64, pattern = "^Rps")
TK64[["percent.RPL"]] <- PercentageFeatureSet(TK64, pattern = "^Rpl")

#filter some rough values
TK64 <- subset(TK64, subset = percent.mt <50)
TK64 <- subset(TK64, subset = percent.RPS <25)

QUALITY <- TK64@meta.data
write.csv(QUALITY, file = "CSV/E11_Rep_QC.csv")

VlnPlot(TK64, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK64 <- SCTransform(TK64, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK64 <- RunPCA(TK64, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK64, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK64 <- RunUMAP(TK64, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK64, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()

####Check####
genes = c("Foxd1","Foxg1","Shh","Fst","Pax6","Sp8","Sp9","Nkx2-1","Rax","Foxa1","Foxa2")
genes <- c("Dlx2", "Olig2", "Nkx2-2", "Olfml3", "Chrd", "Nov", "Car2",
           "Sst", "Fgf18", "Six6", "Dbx1", "Bmp2", "Bmp7", "Vax1", "Vax2", "Emx2", "Nkx2-4", "Sfrp2", "Pitx2")
genes <- c("Six3", "Lef1", "Axin2", "Top2a", "Cdc20", "Ccnb1", "Ccne1")
genes <- c("Sim1", "Otp", "Emx2", "Lhx5", "Lhx1", "Isl1", "Nr5a1", "Wnt8a", "Wnt8b", "Fgf10", "Pax2", "Vax1","Vax2")
genes <- c("Xist","Malat1","Tsix","Rps4x","Ddx3x","Ddx3y")
genes <- c("Avp", "Spp1", "Cst3", "Rorb")
genes <- c("Lhx2", "Lhx9", "Pitx3", "Gbx2", "Tcf7l2", "Gata2", "Gata3", "Hmx2","Hmx3","Nr5a1")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
m.s.genes$mouseGene
m.g2m.genes$mouseGene

p <- FeaturePlot(TK64, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK64 <- CellCycleScoring(TK64, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK64, group.by = "Phase")

#FeaturePlot(TK64, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK64, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK64 <- FindNeighbors(TK64, dims = 1:20)
TK64 <- FindClusters(TK64, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK64, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK64, test.use = "wilcox",
                          logfc.threshold = 0.2,
                          min.pct = 0.01, verbose = T)
write.csv(markers, file = "CSV/E11_Rep_DEG-initial.csv")


####1st Cleanup####
#remove 
#6 Col3a1, meningeal
#8 blood
#11 Col4a1, meningeal

TK64 <- subset(TK64, idents = c("0","1","2","3","4","5",
                                "7","9","10","12","13","14"))
QUALITY <- TK64@meta.data
write.csv(QUALITY, file = "CSV/E11_Rep_QC.csv")
TK64 <- SCTransform(TK64, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
TK64 <- RunPCA(TK64, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK64, ndims = 50, reduction = "pca" )
TK64 <- RunUMAP(TK64, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK64, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()
TK64 <- FindNeighbors(TK64, dims = 1:20)
TK64 <- FindClusters(TK64, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK64, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
markers <- FindAllMarkers(TK64, test.use = "wilcox",
                          logfc.threshold = 0.2,
                          min.pct = 0.01, verbose = T)
write.csv(markers, file = "CSV/E11_Rep_DEG-initial.csv")
#Leave immune cells
#13 blood cells, remove after merging
save(TK64, file = "Robj/E11_Rep.Robj")



#start from here again
#Merge with HyDD_v1 E11
E11 <- merge(x = Hypo2, y = list(TK64))
Idents(E11) <- "orig.ident"
E11[["percent.sex"]] <- PercentageFeatureSet(E11,
                                             features = c("Xist","Malat1","Tsix"))
E11[["percent.mt"]] <- PercentageFeatureSet(E11, pattern = "^mt-")
E11[["percent.RPS"]] <- PercentageFeatureSet(E11, pattern = "^Rps")
E11[["percent.RPL"]] <- PercentageFeatureSet(E11, pattern = "^Rpl")

VlnPlot(E11, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)

table(E11@meta.data$orig.ident)
Idents(E11) <- "Age"
E11<- SCTransform(E11, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA

#Run PCA
E11<- RunPCA(E11, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(E11, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
m.s.genes$mouseGene
m.g2m.genes$mouseGene
E11 <- CellCycleScoring(E11, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)

####harmony####
E11<- RunHarmony(E11, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E11<- RunUMAP(E11, reduction = "harmony", dims = 1:20,
              n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(E11, reduction = "umap", label = F,pt.size = 0.1)

DimPlot(E11, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(E11, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
DimPlot(E11, reduction = "umap", label = F,pt.size = 0.1, group.by = "Cluster" )
DimPlot(E11, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(E11, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
        split.by = "Phase", ncol = 2) + NoLegend()

####Check####
genes = c("Foxd1","Foxg1","Shh","Fst","Pax6","Sp8","Sp9","Nkx2-1","Rax","Foxa1","Foxa2")
genes <- c("Dlx2", "Olig2", "Nkx2-2", "Olfml3", "Chrd", "Nov", "Car2",
           "Sst", "Fgf18", "Six6", "Dbx1", "Bmp2", "Bmp7", "Vax1", "Vax2", "Emx2", "Nkx2-4", "Sfrp2", "Pitx2")
genes <- c("Six3", "Lef1", "Axin2", "Top2a", "Cdc20", "Ccnb1", "Ccne1")
genes <- c("Sim1", "Otp", "Emx2", "Lhx5", "Lhx1", "Isl1", "Nr5a1", "Wnt8a", "Wnt8b", "Fgf10", "Pax2", "Vax1","Vax2")
genes <- c("Xist","Malat1","Tsix","Rps4x","Ddx3x","Ddx3y")
genes <- c("Avp", "Spp1", "Cst3", "Rorb")
genes <- c("Lhx2", "Lhx9", "Pitx3", "Gbx2", "Tcf7l2", "Gata2", "Gata3", "Hmx2","Hmx3")
p <- FeaturePlot(E11, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

####Cluster-identity####
E11 <- FindNeighbors(E11, dims = 1:20, reduction = "harmony")
E11 <- FindClusters(E11, resolution = 2) #SCT_snn_res.0.8
DimPlot(E11, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(E11, file = "Robj/E11_Merge.Robj")
table(E11@active.ident, E11@meta.data$Cluster)

#initial markers
markers <- FindAllMarkers(E11, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/E11-Cluster-DEG-initial.csv")





#genes to plot
genes <- c("Cdkn1c","Gadd45g","Pkib","Dlk1","Hes6","Ass1") #AH Neural Progenitor
genes <- c("D930028M14Rik", "Onecut2", "Sncg", "Snhg11", "Ina", "Cartpt") #Ant ID
genes <- c("Pomc", "Chchd10", "S100a10", "Six6","Rbp1","Sox14") #ARC
genes <- c("Sst","Isl1","Six6","Otp","Dlk1","Cited1") #DMH
genes <- c("Lhx8", "Snhg11","Cadm1","Pcsk1n","Mapt","Gabra2") #ID/TT
genes <- c("Apoe","Lgals1","Igfbp7","Sepp1","Serpinh1","Tyrobp") #Immune cell
genes <- c("Pmch","Lhx9","Chd3os") #LH?
genes <- c("Dlx2","Lhx8","Dlx6os1","Sp9","Arx","Foxg1") #MGE
genes <- c("Lhx1os","Pcp4","Foxb1","Lhx1","Nhlh2","Lhx5") #MMN
genes <- c("Ube2c","Ccnb1","Hmgb2","Spc25","Cks2","Neurog2") #G2M NPC
genes <- c("Pitx2","Gm28050","Nrn1","Uncx","Tmem163","Nhlh2") #Pax7 TT
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

p <- FeaturePlot(E11, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = F, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)



#Rename
#G2M NPC mixed with ZLI
#no Tuberal - ARC/VMH/PMN
#Thalamic progenitor Lhx2 expression is low at this age
#POA not whosn
#SCN not shown
#Weird Sst/Meis2 ergion
#LH not clear

#0 NPC - Tuberal
#1 NPC
#2 NPC Tuberal
#3 Neural progenitor (Ascl1)
#4 G2M NPC?
#5 Neural progenitor (Neurog1/Neurog2)
#6 NPC
#7 SMN/MMN
#8 Prethalamus?
#9 G2M NPC
#10	Foxg1. Junk?
#11 NPC Thalamic?
#12 NPC Thalamic?
#13 NPC Foxg1 
#14 Foxg1. MGE
#15 ZLI or telencephalon roof plate
#16 Foxg1. Thalamic eminence
#17 Foxg1. Thalamic eminence
#18 G2M NPC
#19 PVH & SON
#20 SMN/MMN
#21 G2M NPC
#22 Junk?
#23 G2M NPC
#24 Meis2/Sst
#25 Neural progenitor
#26 Junk?
#27 SMN/MMN or Pax7+ TT
#28 Ant ID/MMN/TT/ID
#29 Immune cells
#30 Immune cells
#31 Immune cells
#32 Foxg1. Junk
#33 Junk
#34 Endothelial remove
#35 Immune cells

E11@meta.data$SCT_snn_res.0.4 <- NULL
E11@meta.data$SCT_snn_res.0.8 <- NULL
E11 <- RenameIdents(E11, "0" = "NPC (Tuberal)", "1" = "NPC", "2" = "NPC (Tuberal)", "3" = "Neural Pro (Ascl1)",
                    "4" = "G2M NPC", "5" = "Neural Pro (Neurog1)", "6" = "NPC", "7" = "SMN/MMN", "8" = "Prethaalmus",
                    "9" = "G2M NPC", "10" = "Junk (Foxg1)", "11" = "NPC (Thalamic)", "12" = "NPC (Thalamic)", "13" = "NPC (Foxg1)",
                    "14" = "MGE", "15" = "Check (Telen RP/ZLI)", "16" = "Thalamic eminence", "17" = "Thalamic eminence",
                    "18" = "G2M NPC", "19" = "PVH/SON", "20" = "SMN/MMN", "21" = "G2M NPC", "22" = "Junk", "23" = "G2M NPC",
                    "24" = "Check (Meis2/Sst)", "25" = "Neural Pro", "26" = "Junk", "27" ="SMN/MMN (TT?)",
                    "28" = "Check (MMN/TT/ID)", "29" = "Immune cells", "30" = "Immune cells", "31" = "Immune cells",
                    "32" = "Junk (Foxg1)", "33" = "Junk", "34" = "Endothelial", "35" = "Immune cells")
Group <- sort(levels(E11@active.ident))
E11 <- AddMetaData(E11, E11@active.ident, "Cluster_Pass1")
E11@active.ident <- factor(x = E11@active.ident, levels = Group)
E11@meta.data$Cluster_Pass1 <- factor(x = E11@meta.data$Cluster_Pass1, levels = Group)
DimPlot(E11, reduction = "umap", label = F, pt.size = 0.5)

save(E11, file = "Robj/E11_Merge.Robj")




