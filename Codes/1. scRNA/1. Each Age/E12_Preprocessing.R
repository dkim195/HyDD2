library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)

####HyDD2####
#Processing E12
#fresh

setwd("D:/")
#Load HyDDV1, subset E12
load(file = "HyDDv1/HyDD_pattern.Robj")
Idents(Hypo2) <- "orig.ident"
Hypo2 <- subset(Hypo2, idents = "E12")
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
write.csv(QUALITY, file = "CSV/E12_HyDD1_QC.csv")
save(Hypo2, file = "Robj/E12_HyDDV1.Robj")

####TK59####
#New v3 10x HyDD2
#TK59 - first replicate
TK59<- Read10X_h5("scRNA/TK59/outs/raw_feature_bc_matrix.h5")
#adding column names to avoid barcode clash
colnames(TK59) = paste0("TK59_", colnames(TK59))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK59 <- CreateSeuratObject(counts = TK59, project = "E12",
                           min.cells = 5, min.features = 1000)
TK59 <- subset(TK59, subset = nCount_RNA > 2000)

###merge####
TK59 <- RenameIdents(TK59, "TK59"= "E12_Rep1")
TK59 <- AddMetaData(TK59, TK59@active.ident, "Age")
head(TK59@meta.data)
Idents(TK59) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK59[["percent.mt"]] <- PercentageFeatureSet(TK59, pattern = "^mt-")
TK59[["percent.RPS"]] <- PercentageFeatureSet(TK59, pattern = "^Rps")
TK59[["percent.RPL"]] <- PercentageFeatureSet(TK59, pattern = "^Rpl")

#filter some rough values
TK59 <- subset(TK59, subset = percent.mt <50)
TK59 <- subset(TK59, subset = percent.RPS <25)

QUALITY <- TK59@meta.data
write.csv(QUALITY, file = "CSV/E12_Rep1_QC.csv")

VlnPlot(TK59, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK59 <- SCTransform(TK59, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK59 <- RunPCA(TK59, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK59, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK59 <- RunUMAP(TK59, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK59, reduction = "umap", label = F,pt.size = 0.5, 
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

p <- FeaturePlot(TK59, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK59 <- CellCycleScoring(TK59, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK59, group.by = "Phase")

#FeaturePlot(TK59, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK59, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK59 <- FindNeighbors(TK59, dims = 1:20)
TK59 <- FindClusters(TK59, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK59, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK59, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E12_Rep1_DEG-initial.csv")


####1st Cleanup####
#remove 
#16 Col4a1
#17 Blood
#keep E18 immune

TK59 <- subset(TK59, idents = c("0","1","2","3","4","5","6","7","8",
                                "9","10","11","12","13","14","15","18"))
QUALITY <- TK59@meta.data
write.csv(QUALITY, file = "CSV/E12_Rep1_QC.csv")
TK59 <- SCTransform(TK59, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
TK59 <- RunPCA(TK59, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK59, ndims = 50, reduction = "pca" )
TK59 <- RunUMAP(TK59, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK59, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()
TK59 <- FindNeighbors(TK59, dims = 1:20)
TK59 <- FindClusters(TK59, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK59, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
markers <- FindAllMarkers(TK59, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E12_Rep1_DEG-initial.csv")
#Leave immune cells
#
save(TK59, file = "Robj/E12_Rep1.Robj")

####TK60####
#New v3 10x HyDD2
#Second repliate
TK60<- Read10X_h5("scRNA/TK60/outs/raw_feature_bc_matrix.h5")
#adding column names to avoid barcode clash
colnames(TK60) = paste0("TK60_", colnames(TK60))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK60 <- CreateSeuratObject(counts = TK60, project = "E12",
                           min.cells = 5, min.features = 1000)
TK60 <- subset(TK60, subset = nCount_RNA > 2000)

###merge####
TK60 <- RenameIdents(TK60, "TK60"= "E12_Rep2")
TK60 <- AddMetaData(TK60, TK60@active.ident, "Age")
head(TK60@meta.data)
Idents(TK60) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK60[["percent.mt"]] <- PercentageFeatureSet(TK60, pattern = "^mt-")
TK60[["percent.RPS"]] <- PercentageFeatureSet(TK60, pattern = "^Rps")
TK60[["percent.RPL"]] <- PercentageFeatureSet(TK60, pattern = "^Rpl")

#filter some rough values
TK60 <- subset(TK60, subset = percent.mt <50)
TK60 <- subset(TK60, subset = percent.RPS <25)

QUALITY <- TK60@meta.data
write.csv(QUALITY, file = "CSV/E12_Rep2_QC.csv")

VlnPlot(TK60, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK60 <- SCTransform(TK60, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK60 <- RunPCA(TK60, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK60, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK60 <- RunUMAP(TK60, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK60, reduction = "umap", label = F,pt.size = 0.5, 
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

p <- FeaturePlot(TK60, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK60 <- CellCycleScoring(TK60, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK60, group.by = "Phase")

#FeaturePlot(TK60, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK60, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK60 <- FindNeighbors(TK60, dims = 1:20)
TK60 <- FindClusters(TK60, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK60, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK60, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E12_Rep2_DEG-initial.csv")

####1st Cleanup####
#remove 
#19 blood

TK60 <- subset(TK60, idents = c("0","1","2","3","4","5","6","7","8",
                                "9","10","11","12","13","14","15","16","17","18"))
QUALITY <- TK60@meta.data
write.csv(QUALITY, file = "CSV/E12_Rep2_QC.csv")
TK60 <- SCTransform(TK60, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
TK60 <- RunPCA(TK60, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK60, ndims = 50, reduction = "pca" )
TK60 <- RunUMAP(TK60, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK60, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()
TK60 <- FindNeighbors(TK60, dims = 1:20)
TK60 <- FindClusters(TK60, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK60, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
markers <- FindAllMarkers(TK60, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E12_Rep2_DEG-initial.csv")
#Leave immune cells
#13 blood cells, remove after merging
save(TK60, file = "Robj/E12_Rep2.Robj")



load(file = "Robj/E12_HyDDV1.Robj")
#start from here again
#Merge with HyDD_v1 E12
E12 <- merge(x = Hypo2, y = list(TK59, TK60))
head(E12@meta.data)
E12@meta.data$percent.mt <- NULL
E12@meta.data$percent.RPS <- NULL
E12@meta.data$percent.RPL <- NULL
E12@meta.data$S.Score <- NULL
E12@meta.data$G2M.Score <- NULL
E12@meta.data$Phase <- NULL
E12@meta.data$old.ident <- NULL
E12@meta.data$SCT_snn_res.0.4 <- NULL
E12@meta.data$seurat_clusters <- NULL
save(E12, file = "Robj/E12_Merge.Robj")
Idents(E12) <- "orig.ident"
E12[["percent.sex"]] <- PercentageFeatureSet(E12,
                                             features = c("Xist","Malat1","Tsix"))
E12[["percent.mt"]] <- PercentageFeatureSet(E12, pattern = "^mt-")
E12[["percent.RPS"]] <- PercentageFeatureSet(E12, pattern = "^Rps")
E12[["percent.RPL"]] <- PercentageFeatureSet(E12, pattern = "^Rpl")


VlnPlot(E12, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                          "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)

table(E12@meta.data$orig.ident)
Idents(E12) <- "Age"
E12<- SCTransform(E12, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA

#Run PCA
E12<- RunPCA(E12, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(E12, ndims = 50, reduction = "pca" )
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
E12 <- CellCycleScoring(E12, s.features = m.s.genes$mouseGene,
                        g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)

####harmony####
E12<- RunHarmony(E12, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E12<- RunUMAP(E12, reduction = "harmony", dims = 1:20,
              n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(E12, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(E12, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
DimPlot(E12, reduction = "umap", label = F,pt.size = 0.1, group.by = "Cluster" )
DimPlot(E12, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(E12, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
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
p <- FeaturePlot(E12, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

save(E12, file = "Robj/E12_Merge.Robj")


####Cluster-identity####
E12 <- FindNeighbors(E12, dims = 1:20, reduction = "harmony")
E12 <- FindClusters(E12, resolution = 2) #SCT_snn_res.0.8
DimPlot(E12, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
table(E12@active.ident, E12@meta.data$Cluster)

#initial markers
markers <- FindAllMarkers(E12, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/E12-Cluster-DEG-initial.csv")





#genes to plot
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

p <- FeaturePlot(E12, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = F, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)



#Rename
#G2M NPC mixed with ZLI
#Tuberal shown
#Thalamic progenitor Lhx2 high
#POA not shown but Foxg1/PVH/SON
#SCN not shown but Foxg1/PVH/SON
#LH not clear
#Prethalamus high development

#0 SMN
#1 Prethalamus
#2 G2M NPC
#3 Junk?
#4 Prethalamus
#5 PVH/SON
#6 Check - Junk?
#7 SMN
#8 MMN
#9 NPC
#10 NPC
#11 PVH/SON
#12 NPC - ASCL1 / GABAergic?
#13 G2M NPC or NPC?
#14 NPC -Neurog2 / SMN/MMN pro
#15 Junk
#16 Tuberal (PMN)
#17 Prethalamus (has some ID)
#18 Foxg1/PVH/SON
#19 Check - Junk?
#20 Thalamic eminence
#21 G2M NPC or NPC?
#22 Tuberal (ARC/VMH?)
#23 Tuberal (ARC/VMH?)
#24 Tuberal (PMN - I see LH)
#25 Prethalamus
#26 MMN
#27 SMN/MMN progenitor
#28 Prethalamus
#29 G2M NPC or NPC?
#30 G2M NPC or NPC?
#31 Foxg1 (NPC/Thalamic eminence?)
#32 G2M NPC or NPC?
#33 G2M NPC
#34 G2M NPC or NPC?
#35 Prethalamus
#36 Prethalamus (ID)
#37 SMN
#38 G2M NPC or NPC
#39 Immune cells
#40 Rim domain
#41 NPC (Neurog2)
#42 Immune cells



E12 <- RenameIdents(E12, "0" = "SMN", "1" = "Prethalamus", "2" = "G2M NPC/NPC", "3" = "Junk", "4" = "Prethalamus",
                    "5" = "PVH/SON", "6" = "Check (Junk?)", "7" = "SMN", "8" = "MMN", "9" = "NPC", "10" = "NPC",
                    "11" = "PVH/SON", "12" = "Neural Pro (Ascl1)", "13" = "G2M NPC/NPC", "14" = "Neural Pro (Neurog2)",
                    "15" = "Junk", "16" = "Tuberal (PMN)", "17" = "Prethalamus (ID)", "18" = "Foxg1/PVH/SON",
                    "19" = "Check (Junk?)", "20" = "Thalamic eminence", "21" = "G2M NPC/NPC", "22" = "Tuberal (ARC/VMH)",
                    "23" ="Tuberal (ARC/VMH)", "24" = "Tuberal (PMN-LH)", "25" = "Prethalamus", "26" = "MMN",
                    "27" = "Neural Pro (SMN/MMN)", "28" = "Prethalamus", "29" = "G2M NPC/NPC", "30" = "G2M NPC/NPC",
                    "31" = "Foxg1 (NPC/Thalamic eminence?)", "32" = "G2M NPC/NPC", "33" = "G2M NPC/NPC", "34" = "G2M NPC/NPC",
                    "35" = "Prethalamus", "36" = "Prethalamus (ID)", "37" = "SMN", 
                    "38" = "G2M NPC/NPC", "39" = "Immune cells",
                    "40" = "Rim domain", "41" = "Neural Pro (Neurog2)", "42" = "Immune cells")
Group <- sort(levels(E12@active.ident))
E12 <- AddMetaData(E12, E12@active.ident, "Cluster_Pass1")
E12@active.ident <- factor(x = E12@active.ident, levels = Group)
E12@meta.data$Cluster_Pass1 <- factor(x = E12@meta.data$Cluster_Pass1, levels = Group)
DimPlot(E12, reduction = "umap", label = F, pt.size = 0.5)

save(E12, file = "Robj/E12_Merge.Robj")


