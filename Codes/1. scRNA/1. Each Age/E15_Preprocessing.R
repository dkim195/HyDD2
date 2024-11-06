library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("D:/")
####HyDD2####
#E15.5 
#Fresh
####TK67####
#New v3.1 10x HyDD2
#TK67 - first replicate
TK67<- Read10X_h5("scRNA/TK67/outs/raw_feature_bc_matrix.h5")
#adding column names to avoid barcode clash
colnames(TK67) = paste0("TK67_", colnames(TK67))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK67 <- CreateSeuratObject(counts = TK67, project = "E15",
                           min.cells = 5, min.features = 1000)
TK67 <- subset(TK67, subset = nCount_RNA > 2000)

###merge####
TK67 <- RenameIdents(TK67, "TK67"= "E15_Rep1")
TK67 <- AddMetaData(TK67, TK67@active.ident, "Age")
head(TK67@meta.data)
Idents(TK67) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK67[["percent.mt"]] <- PercentageFeatureSet(TK67, pattern = "^mt-")
TK67[["percent.RPS"]] <- PercentageFeatureSet(TK67, pattern = "^Rps")
TK67[["percent.RPL"]] <- PercentageFeatureSet(TK67, pattern = "^Rpl")

#filter some rough values
TK67 <- subset(TK67, subset = percent.mt <50)
TK67 <- subset(TK67, subset = percent.RPS <25)

QUALITY <- TK67@meta.data
write.csv(QUALITY, file = "CSV/E15_Rep1_QC.csv")

VlnPlot(TK67, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK67 <- SCTransform(TK67, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK67 <- RunPCA(TK67, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK67, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK67 <- RunUMAP(TK67, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK67, reduction = "umap", label = F,pt.size = 0.5, 
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
genes <- c("Foxd1","Lhx8","Shh","Fst","Pax6","Sp8","Sp9","Nkx2-1","Rax","Foxa1","Foxa2","Lhx9","Sim1","Otp","Nr5a1",
           "Agrp","Npy","Gal","Slc32a1","Slc17a6","Pomc","Cartpt","Hcrt","Nts","Sst","Ghrh","Gnrh1","Lepr",
           "Avp","Crh","Oxt","Th","Vip","Trh","Rfrp")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
m.s.genes$mouseGene
m.g2m.genes$mouseGene

p <- FeaturePlot(TK67, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK67 <- CellCycleScoring(TK67, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK67, group.by = "Phase")

#FeaturePlot(TK67, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK67, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK67 <- FindNeighbors(TK67, dims = 1:20)
TK67 <- FindClusters(TK67, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK67, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK67, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E15_Rep1_DEG-initial.csv")


####1st Cleanup####
#remove 
#11 Blood

TK67 <- subset(TK67, idents = c("0","1","2","3","4","5","6","7","8",
                                "9","10","12","13","14","15"))

QUALITY <- TK67@meta.data
write.csv(QUALITY, file = "CSV/E15_Rep1_QC.csv")
TK67 <- SCTransform(TK67, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
TK67 <- RunPCA(TK67, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK67, ndims = 50, reduction = "pca" )
TK67 <- RunUMAP(TK67, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK67, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()
TK67 <- FindNeighbors(TK67, dims = 1:20)
TK67 <- FindClusters(TK67, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK67, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
markers <- FindAllMarkers(TK67, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E15_Rep1_DEG-initial.csv")
#Leave immune cells
#
save(TK67, file = "Robj/E15_Rep1.Robj")

####TK68####
#New v3.1 10x HyDD2
#Second repliate
TK68<- Read10X_h5("scRNA/TK68/outs/raw_feature_bc_matrix.h5")
#adding column names to avoid barcode clash
colnames(TK68) = paste0("TK68_", colnames(TK68))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK68 <- CreateSeuratObject(counts = TK68, project = "E15",
                           min.cells = 5, min.features = 1000)
TK68 <- subset(TK68, subset = nCount_RNA > 2000)

###merge####
TK68 <- RenameIdents(TK68, "TK68"= "E15_Rep2")
TK68 <- AddMetaData(TK68, TK68@active.ident, "Age")
head(TK68@meta.data)
Idents(TK68) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK68[["percent.mt"]] <- PercentageFeatureSet(TK68, pattern = "^mt-")
TK68[["percent.RPS"]] <- PercentageFeatureSet(TK68, pattern = "^Rps")
TK68[["percent.RPL"]] <- PercentageFeatureSet(TK68, pattern = "^Rpl")
#mito.genes <- grep(pattern = "^mt-", x = rownames(x = TK68@assays$SCT@data), value = TRUE)
#percent.mito <- colSums(expm1(TK68@assays$SCT@data[mito.genes, ]))/colSums(expm1(TK68@assays$SCT@data))
#TK68 <- AddMetaData(TK68, percent.mito, "percent.mito")

#filter some rough values
TK68 <- subset(TK68, subset = percent.mt <50)
TK68 <- subset(TK68, subset = percent.RPS <25)

QUALITY <- TK68@meta.data
write.csv(QUALITY, file = "CSV/E15_Rep2_QC.csv")

VlnPlot(TK68, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK68 <- SCTransform(TK68, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK68 <- RunPCA(TK68, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK68, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK68 <- RunUMAP(TK68, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK68, reduction = "umap", label = F,pt.size = 0.5, 
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

p <- FeaturePlot(TK68, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK68 <- CellCycleScoring(TK68, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK68, group.by = "Phase")

#FeaturePlot(TK68, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK68, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK68 <- FindNeighbors(TK68, dims = 1:20)
TK68 <- FindClusters(TK68, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK68, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK68, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E15_Rep2_DEG-initial.csv")

####1st Cleanup####
#remove 
#13 haemoglobin

TK68 <- subset(TK68, idents = c("0","1","2","3","4","5","6","7","8",
                                "9","10","11","12","14","15"))
QUALITY <- TK68@meta.data
write.csv(QUALITY, file = "CSV/E15_Rep2_QC.csv")
TK68 <- SCTransform(TK68, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
TK68 <- RunPCA(TK68, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK68, ndims = 50, reduction = "pca" )
TK68 <- RunUMAP(TK68, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK68, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()
TK68 <- FindNeighbors(TK68, dims = 1:20)
TK68 <- FindClusters(TK68, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK68, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
markers <- FindAllMarkers(TK68, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E15_Rep2_DEG-initial.csv")
#Leave immune cells
save(TK68, file = "Robj/E15_Rep2.Robj")


#start from here again
#Merge with HyDD_v1 E15
E15 <- merge(x = TK67, y = list(TK68))
head(E15@meta.data)
E15@meta.data$percent.mt <- NULL
E15@meta.data$percent.RPS <- NULL
E15@meta.data$percent.RPL <- NULL
E15@meta.data$S.Score <- NULL
E15@meta.data$G2M.Score <- NULL
E15@meta.data$Phase <- NULL
E15@meta.data$old.ident <- NULL
E15@meta.data$SCT_snn_res.0.4 <- NULL
E15@meta.data$seurat_clusters <- NULL
save(E15, file = "Robj/E15_Merge.Robj")
load(file = "Robj/E15_Merge.Robj")
Idents(E15) <- "orig.ident"
E15[["percent.sex"]] <- PercentageFeatureSet(E15,
                                             features = c("Xist","Malat1","Tsix"))
E15[["percent.mt"]] <- PercentageFeatureSet(E15, pattern = "^mt-")
E15[["percent.RPS"]] <- PercentageFeatureSet(E15, pattern = "^Rps")
E15[["percent.RPL"]] <- PercentageFeatureSet(E15, pattern = "^Rpl")


VlnPlot(E15, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                          "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)

table(E15@meta.data$orig.ident)
Idents(E15) <- "Age"
E15<- SCTransform(E15, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA

#Run PCA
E15<- RunPCA(E15, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(E15, ndims = 50, reduction = "pca" )
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
E15 <- CellCycleScoring(E15, s.features = m.s.genes$mouseGene,
                        g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)

####harmony####
E15<- RunHarmony(E15, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E15<- RunUMAP(E15, reduction = "harmony", dims = 1:20,
              n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(E15, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(E15, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
DimPlot(E15, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(E15, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
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
genes <- c("Foxd1","Lhx8","Shh","Fst","Pax6","Sp8","Sp9","Nkx2-1","Rax","Foxa1","Foxa2","Lhx9","Sim1","Otp",
           "Agrp","Npy","Gal","Slc32a1","Slc17a6","Pomc","Cartpt","Hcrt","Nts","Sst","Ghrh","Gnrh1","Lepr",
           "Avp","Crh","Oxt","Th","Vip","Trh","Rfrp","Lhx6","Foxg1","Nkx2-2","Pvalb","Tac1")
p <- FeaturePlot(E15, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

save(E15, file = "Robj/E15_Merge.Robj")

####Cluster-identity####
E15 <- FindNeighbors(E15, dims = 1:20, reduction = "harmony")
E15 <- FindClusters(E15, resolution = 2) #SCT_snn_res.0.8
DimPlot(E15, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(E15, file = "Robj/E15_Merge.Robj")
table(E15@active.ident, E15@meta.data$Cluster)

#initial markers
markers <- FindAllMarkers(E15, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/E15-Cluster-DEG-initial.csv")





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

p <- FeaturePlot(E15, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = F, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)



#Rename
#E15
#- Pattern still maintain
#- Better separation of unidentified/hidden areas?
#  = where is POA/SCN/DMH

#0 Prethalamus? ID?
#1 PMN
#2 MMN
#3 ARC?
#4 neurogenic? gli genic?
#5 PVH/SON?
#6 SMN
#7 Prethalamus? ID?
#8 VMH
#9 MMN
#10 SMN
#11 VMH
#12 SMN
#13 Foxg1 Neural progenitor?
#14 PVH/SON?
#15 MMN CHECK
#16 Foxg1 Neural progenitor?
#17 Neural progenitor? GLIAL?
#18 POMC
#19 Junk?
#20 Prethal (GABAergic?)
#21 SMN/MMN progenitor
#22 SMN
#23 MMN (CCK)
#24 Junk
#25 ARC?
#26 Prethalamus (oligo?)
#27 POA/SCN
#28 Neural progenitor? GLIAL?
#29 Foxg1 Neural progenitor? GLIAL?
#30 SMN (glial?)
#31 ID?
#32 ID?
#33 Foxg1 Neural progenitor? GLIAL?
#34 ARC (SST/NPY)
#35 ARC (glial?)
#36 Foxg1 Junk?
#37 PVN/SON?
#38 Neural progenitor? GLIAL?
#39 Foxg1 ascl1+ neurogenic? gli genic?
#40 endothelial junk
#41 Foxg1 endothelial junk
#42 ARC?
#43 Pmch neurons
#44 ID?
#45 Thalamic?
#46 Immune
#47 Hcrt/Npvf
#48 Foxg1 Junk


E15 <- RenameIdents(E15, "0" = "PreThal/ID", "1" = "PMM", "2" = "MMN", "3" = "ARC?", "4" = "Glio Pro?",
                    "5" = "PVH/SON?", "6" = "SMN", "7" = "PreThal/ID", "8" = "VMH", "9" = "MMN",
                    "10" = "SMN", "11" = "VMH", "12" = "SMN", "13" = "Foxg1", "14" = "PVH/SON?",
                    "15" = "Check (MMN?)", "16" = "Foxg1", "17" = "Glio Pro?", "18" = "Pomc", "19" = "Junk",
                    "20" = "PreThal/ID", "21" = "Neural Pro (SMN/MMN)", "22" = "SMN", "23" = "MMN (Cck)",
                    "24" = "Junk", "25" = "Arc?", "26" = "PreThal/ID (Glia?)", "27" = "Check (POA/SCN)",
                    "28" = "Glio Pro?", "29" = "Foxg1", "30" = "SMN (Glia?)", "31" = "Check (ID?)",
                    "32" = "Check (ID?)", "33" = "Foxg1", "34" = "ARC (Sst/Npy)", "35" = "Arc (Glial?)",
                    "36" = "Foxg1", "37" = "PVH/SON?", "38" = "Glio Pro?", "39" = "Foxg1 (Ascl1+)", 
                    "40"  = "Endothelial", "41" = "Foxg1", "42" = "Arc?", "43" = "Pmch", "44" = "Check (ID?)",
                    "45" = "Check (Thalamic?)", "46" = "Immune", "47" = "Hcrt/Npvf", "48" = "Junk")
Group <- sort(levels(E15@active.ident))
E15 <- AddMetaData(E15, E15@active.ident, "Cluster_Pass1")
E15@active.ident <- factor(x = E15@active.ident, levels = Group)
E15@meta.data$Cluster_Pass1 <- factor(x = E15@meta.data$Cluster_Pass1, levels = Group)
DimPlot(E15, reduction = "umap", label = F, pt.size = 0.5)

save(E15, file = "Robj/E15_Merge.Robj")






