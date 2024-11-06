library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)

####HyDD2####
#Processing E13
#fresh

setwd("D:/")
#Load HyDDV1, subset E13
load(file = "HyDDv1/HyDD_pattern.Robj")
Idents(Hypo2) <- "orig.ident"
Hypo2 <- subset(Hypo2, idents = "E13")
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
write.csv(QUALITY, file = "CSV/E13_HyDD1_QC.csv")
save(Hypo2, file = "Robj/E13_HyDDV1.Robj")

####TKX####
#New v3 10x HyDD2
#TKX - first replicate
TKX<- Read10X_h5("scRNA/TKX/outs/raw_feature_bc_matrix.h5")
#adding column names to avoid barcode clash
colnames(TKX) = paste0("TKX_", colnames(TKX))
TKX <- CreateSeuratObject(counts = TKX, project = "E13",
                           min.cells = 5, min.features = 1000)
TKX <- subset(TKX, subset = nCount_RNA > 2000) 

###merge####
TKX <- RenameIdents(TKX, "TKX"= "E13_Rep1")
TKX <- AddMetaData(TKX, TKX@active.ident, "Age")
head(TKX@meta.data)
Idents(TKX) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TKX[["percent.mt"]] <- PercentageFeatureSet(TKX, pattern = "^mt-")
TKX[["percent.RPS"]] <- PercentageFeatureSet(TKX, pattern = "^Rps")
TKX[["percent.RPL"]] <- PercentageFeatureSet(TKX, pattern = "^Rpl")

#filter some rough values
TKX <- subset(TKX, subset = percent.mt <50)
TKX <- subset(TKX, subset = percent.RPS <25)

QUALITY <- TKX@meta.data
write.csv(QUALITY, file = "CSV/E13_Rep1_QC.csv")

VlnPlot(TKX, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TKX <- SCTransform(TKX, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TKX <- RunPCA(TKX, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TKX, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TKX <- RunUMAP(TKX, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TKX, reduction = "umap", label = F,pt.size = 0.5, 
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

p <- FeaturePlot(TKX, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TKX <- CellCycleScoring(TKX, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TKX, group.by = "Phase")

#FeaturePlot(TKX, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TKX, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TKX <- FindNeighbors(TKX, dims = 1:20)
TKX <- FindClusters(TKX, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TKX, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TKX, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E13_Rep1_DEG-initial.csv")


####1st Cleanup####
#remove 
#12 blood
TKX <- subset(TKX, idents = c("0","1","2","3","4","5","6","7","8",
                              "9","10","11","13","14","15"))
QUALITY <- TKX@meta.data
write.csv(QUALITY, file = "CSV/E13_Rep1_QC.csv")
TKX <- SCTransform(TKX, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
TKX <- RunPCA(TKX, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TKX, ndims = 50, reduction = "pca" )
TKX <- RunUMAP(TKX, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TKX, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()
TKX <- FindNeighbors(TKX, dims = 1:20)
TKX <- FindClusters(TKX, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TKX, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
markers <- FindAllMarkers(TKX, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E13_Rep1_DEG-initial.csv")
save(TKX, file = "Robj/E13_Rep1.Robj")



####TK100####
#New v3 10x HyDD2
#Second repliate
TK100<- Read10X_h5("scRNA/TK100/outs/raw_feature_bc_matrix.h5")
#adding column names to avoid barcode clash
colnames(TK100) = paste0("TK100_", colnames(TK100))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK100 <- CreateSeuratObject(counts = TK100, project = "E13",
                           min.cells = 5, min.features = 1000)

TK100 <- subset(TK100, subset = nCount_RNA > 2000) 

###merge####
TK100 <- RenameIdents(TK100, "TK100"= "E13_Rep2")
TK100 <- AddMetaData(TK100, TK100@active.ident, "Age")
head(TK100@meta.data)
Idents(TK100) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK100[["percent.mt"]] <- PercentageFeatureSet(TK100, pattern = "^mt-")
TK100[["percent.RPS"]] <- PercentageFeatureSet(TK100, pattern = "^Rps")
TK100[["percent.RPL"]] <- PercentageFeatureSet(TK100, pattern = "^Rpl")

#filter some rough values
TK100 <- subset(TK100, subset = percent.mt <50)
TK100 <- subset(TK100, subset = percent.RPS <25)

QUALITY <- TK100@meta.data
write.csv(QUALITY, file = "CSV/E13_Rep2_QC.csv")

VlnPlot(TK100, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK100 <- SCTransform(TK100, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK100 <- RunPCA(TK100, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK100, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK100 <- RunUMAP(TK100, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK100, reduction = "umap", label = F,pt.size = 0.5, 
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

p <- FeaturePlot(TK100, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK100 <- CellCycleScoring(TK100, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK100, group.by = "Phase")

#FeaturePlot(TK100, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK100, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK100 <- FindNeighbors(TK100, dims = 1:20)
TK100 <- FindClusters(TK100, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK100, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK100, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E13_Rep2_DEG-initial.csv")

####1st Cleanup####
#remove 
#14 blood cell

TK100 <- subset(TK100, idents = c("0","1","2","3","4","5","6","7","8","9",
                                  "10","11","12","13","15","16","17","18"))
QUALITY <- TK100@meta.data
write.csv(QUALITY, file = "CSV/E13_Rep2_QC.csv")
TK100 <- SCTransform(TK100, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
TK100 <- RunPCA(TK100, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK100, ndims = 50, reduction = "pca" )
TK100 <- RunUMAP(TK100, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK100, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()
TK100 <- FindNeighbors(TK100, dims = 1:20)
TK100 <- FindClusters(TK100, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK100, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
markers <- FindAllMarkers(TK100, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E13_Rep2_DEG-initial.csv")
#Leave immune cells
save(TK100, file = "Robj/E13_Rep2.Robj")



load(file = "Robj/E13_HyDDV1.Robj")
#start from here again
#Merge with HyDD_v1 E13
E13 <- merge(x = Hypo2, y = list(TKX, TK100))
head(E13@meta.data)
E13@meta.data$percent.mt <- NULL
E13@meta.data$percent.RPS <- NULL
E13@meta.data$percent.RPL <- NULL
E13@meta.data$S.Score <- NULL
E13@meta.data$G2M.Score <- NULL
E13@meta.data$Phase <- NULL
E13@meta.data$old.ident <- NULL
E13@meta.data$SCT_snn_res.0.4 <- NULL
E13@meta.data$seurat_clusters <- NULL
save(E13, file = "Robj/E13_Merge.Robj")
Idents(E13) <- "orig.ident"
E13[["percent.sex"]] <- PercentageFeatureSet(E13,
                                             features = c("Xist","Malat1","Tsix"))
E13[["percent.mt"]] <- PercentageFeatureSet(E13, pattern = "^mt-")
E13[["percent.RPS"]] <- PercentageFeatureSet(E13, pattern = "^Rps")
E13[["percent.RPL"]] <- PercentageFeatureSet(E13, pattern = "^Rpl")


VlnPlot(E13, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                          "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)

table(E13@meta.data$orig.ident)
Idents(E13) <- "Age"
E13<- SCTransform(E13, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA

#Run PCA
E13<- RunPCA(E13, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(E13, ndims = 50, reduction = "pca" )
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
E13 <- CellCycleScoring(E13, s.features = m.s.genes$mouseGene,
                        g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)

####harmony####
E13<- RunHarmony(E13, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E13<- RunUMAP(E13, reduction = "harmony", dims = 1:20,
              n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(E13, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(E13, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
DimPlot(E13, reduction = "umap", label = F,pt.size = 0.1, group.by = "Cluster" )
DimPlot(E13, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(E13, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
        split.by = "Phase", ncol = 2) + NoLegend()

####Check####
genes = c("Foxd1","Foxg1","Shh","Fst","Pax6","Sp8","Sp9","Nkx2-1","Rax","Foxa1","Foxa2")
genes <- c("Dlx2", "Olig2", "Nkx2-2", "Olfml3", "Chrd", "Nov", "Car2",
           "Sst", "Fgf18", "Six6", "Dbx1", "Bmp2", "Bmp7", "Vax1", "Vax2", "Emx2", "Nkx2-4", "Sfrp2", "Pitx2")
genes <- c("Six3", "Lef1", "Axin2", "Top2a", "Cdc20", "Ccnb1", "Ccne1")
genes <- c("Sim1", "Otp", "Emx2", "Lhx5", "Lhx1", "Isl1", "Nr5a1", "Wnt8a", "Wnt8b", "Fgf10", "Pax2", "Vax1","Vax2")
genes <- c("Agrp","Pomc","Pmch","Gal","Cartpt","Npy")
genes <- c("Xist","Malat1","Tsix","Rps4x","Ddx3x","Ddx3y")
genes <- c("Avp", "Spp1", "Cst3", "Rorb")
genes <- c("Lhx2", "Lhx9", "Pitx3", "Gbx2", "Tcf7l2", "Gata2", "Gata3", "Hmx2","Hmx3")
p <- FeaturePlot(E13, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

save(E13, file = "Robj/E13_Merge.Robj")





####Cluster-identity####
E13 <- FindNeighbors(E13, dims = 1:20, reduction = "harmony")
E13 <- FindClusters(E13, resolution = 2) #SCT_snn_res.0.8
DimPlot(E13, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
table(E13@active.ident)

#initial markers
markers <- FindAllMarkers(E13, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/E13-Cluster-DEG-initial.csv")


#genes to plot
genes <- c("Cdkn1c","Gadd45g","Pkib","Dlk1","Hes6","Ass1") #AH Neural Progenitor
genes <- c("D930028M14Rik", "Onecut2", "Sncg", "Snhg11", "Ina", "Cartpt") #Ant ID
genes <- c("Pomc", "Chchd10", "S100a10", "Six6","Rbp1","Sox14") #ARC
genes <- c("Sst","Isl1","Six6","Otp","Dlk1","Cited1") #DMH
genes <- c("Lhx8", "Snhg11","Cadm1","Pcsk1n","Mapt","Gabra2","Lhx6") #ID/TT
genes <- c("Apoe","Lgals1","Igfbp7","Sepp1","Serpinh1","Tyrobp") #Immune cell
genes <- c("Pmch","Lhx9","Chd3os","Hcrt") #LH?
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

p <- FeaturePlot(E13, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = F, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)



#Rename


#G2M NPC mixed with ZLI

#AntID ID clear.
#Foxg1/PVH/SON still there
#Otp+ MMN
#scn HIDDEN IN ANT-/D and pvh/son?


#0 Junk (Thalamic)
#1 Prethalamus
#2 Prethalamus/Foxg1
#3 MMN
#4 PMN
#5 NPC
#6 SMN
#7 NPC
#8 SMN
#9 NEU PRO ASCL1
#10 NEU PRO NEUROG2
#11 Prethalamus
#12 Tuberal (ARC/VMH)
#13 NPC
#14 PHV/SON
#15 NPC
#16 NPC
#17 NPC
#18 Junk
#19 SMN/MMN progenitor
#20 PHV/SON
#21 Tuberal (ARC/VMH)
#22 Thalamic eminence
#23 NEU PRO NEUROG2 (SMN/MMN)
#24 PVH & SON
#25 Foxg1/PVH & SON
#26 JUNK
#27 NEU PRO NEUROG2
#28 NPC
#29 Junk
#30 NPC
#31 Foxg1
#32 Check (MMN?)
#33 ID
#34 Prethalamus (Ant ID?)
#35 Foxg1 NPC
#36 Check (Junk?)
#37 NPC
#38 NPC
#39 NPC
#40 Immune cells
#41 Check (Pax7 TT)
#42 Ant ID
#43 Immune cells
#44 NPC
#45 Tuberal (PMN)
#46 Tuberal (PMN-LH)
#47 jUNK
#48 NEU PRO NEUROG2
#49 NEU PRO NEUROG2
#50 iMMUE CELL



E13 <- RenameIdents(E13, "0" = "Junk (Thalamic)", "1" = "Prethalamus", "2" = "Check (MGE)", "3" = "MMN",
                    "4" = "PMN", "5" = "NPC", "6" = "SMN", "7" = "NPC", "8" = "SMN", "9" = "Neural Pro (Ascl1)",
                    "10" = "Neural Pro (Neurog2)", "11" = "Prethalamus", "12" = "Tuberal (ARC/VMH)", "13" = "NPC",
                    "14" = "PVH/SON", "15" = "NPC", "16" = "NPC", "17" = "NPC", "18" = "Junk",
                    "19" = "Neural Pro (SMN/MMN)",
                    "20" = "PVH/SON", "21" = "Tuberal (ARC/VMH)", "22" = "Thalamic eminence",
                    "23" = "Neural Pro (Neurog2/SMN/MMN)",
                    "24" = "PVH/SON", "25" = "Foxg1/PVH/SON", "26" = "Junk", "27" = "Neural Pro (Neurog2)",
                    "28" = "NPC", "29" = "Junk", "30" = "NPC", "31" = "Foxg1/PVH/SON", "32" = "Check (MMN)",
                    "33" = "ID?", "34" = "Prethalamus (Ant ID?)", "35" = "NPC (Foxg1)", "36" = "Check (Junk?)",
                    "37" = "NPC", "38" = "NPC", "39" = "NPC", "40" = "Immune cells", "41" = "Check (Pax7 TT?)",
                    "42" = "Ant ID?", "43" = "Immune cells", "44" = "NPC", "45" = "Tuberal (PMN)",
                    "46" = "Tuberal (PMN-LH)", "47" = "Junk", "48" = "Neural Pro (Neurog2)",
                    "49" = "Neural Pro (Neurog2)", "50" = "Immune cells")

Group <- sort(levels(E13@active.ident))
E13 <- AddMetaData(E13, E13@active.ident, "Cluster_Pass1")
E13@active.ident <- factor(x = E13@active.ident, levels = Group)
E13@meta.data$Cluster_Pass1 <- factor(x = E13@meta.data$Cluster_Pass1, levels = Group)
DimPlot(E13, reduction = "umap", label = F, pt.size = 0.5)

save(E13, file = "Robj/E13_Merge.Robj")
