library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("D:/")
####HyDD2####
#E14.5 
#Fresh
####TK65####
#New v3.1 10x HyDD2
#TK65 - first replicate
TK65<- Read10X_h5("scRNA/TK65/outs/raw_feature_bc_matrix.h5")
#adding column names to avoid barcode clash
colnames(TK65) = paste0("TK65_", colnames(TK65))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK65 <- CreateSeuratObject(counts = TK65, project = "E14",
                           min.cells = 5, min.features = 1000)
TK65 <- subset(TK65, subset = nCount_RNA > 2000)

###merge####
TK65 <- RenameIdents(TK65, "TK65"= "E14_Rep1")
TK65 <- AddMetaData(TK65, TK65@active.ident, "Age")
head(TK65@meta.data)
Idents(TK65) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK65[["percent.mt"]] <- PercentageFeatureSet(TK65, pattern = "^mt-")
TK65[["percent.RPS"]] <- PercentageFeatureSet(TK65, pattern = "^Rps")
TK65[["percent.RPL"]] <- PercentageFeatureSet(TK65, pattern = "^Rpl")

#filter some rough values
TK65 <- subset(TK65, subset = percent.mt <50)
TK65 <- subset(TK65, subset = percent.RPS <25)

QUALITY <- TK65@meta.data
write.csv(QUALITY, file = "CSV/E14_Rep1_QC.csv")

VlnPlot(TK65, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK65 <- SCTransform(TK65, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK65 <- RunPCA(TK65, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK65, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK65 <- RunUMAP(TK65, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK65, reduction = "umap", label = F,pt.size = 0.5, 
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

p <- FeaturePlot(TK65, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK65 <- CellCycleScoring(TK65, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK65, group.by = "Phase")

#FeaturePlot(TK65, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK65, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK65 <- FindNeighbors(TK65, dims = 1:20)
TK65 <- FindClusters(TK65, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK65, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK65, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E14_Rep1_DEG-initial.csv")


####1st Cleanup####
#remove 
#14 blood
TK65 <- subset(TK65, idents = c("0","1","2","3","4","5","6","7","8",
                                "9","10","11","12","13","15","16"))
QUALITY <- TK65@meta.data
write.csv(QUALITY, file = "CSV/E14_Rep1_QC.csv")
TK65 <- SCTransform(TK65, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
TK65 <- RunPCA(TK65, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK65, ndims = 50, reduction = "pca" )
TK65 <- RunUMAP(TK65, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK65, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()
TK65 <- FindNeighbors(TK65, dims = 1:20)
TK65 <- FindClusters(TK65, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK65, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
markers <- FindAllMarkers(TK65, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E14_Rep1_DEG-initial.csv")
#Leave immune cells
#
save(TK65, file = "Robj/E14_Rep1.Robj")

####TK66####
#New v3.1 10x HyDD2
#Second repliate
TK66<- Read10X_h5("scRNA/TK66/outs/raw_feature_bc_matrix.h5") 
#adding column names to avoid barcode clash
colnames(TK66) = paste0("TK66_", colnames(TK66))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK66 <- CreateSeuratObject(counts = TK66, project = "E14",
                           min.cells = 5, min.features = 1000) 
TK66 <- subset(TK66, subset = nCount_RNA > 2000)

###merge####
TK66 <- RenameIdents(TK66, "TK66"= "E14_Rep2")
TK66 <- AddMetaData(TK66, TK66@active.ident, "Age")
head(TK66@meta.data)
Idents(TK66) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK66[["percent.mt"]] <- PercentageFeatureSet(TK66, pattern = "^mt-")
TK66[["percent.RPS"]] <- PercentageFeatureSet(TK66, pattern = "^Rps")
TK66[["percent.RPL"]] <- PercentageFeatureSet(TK66, pattern = "^Rpl")

#filter some rough values
TK66 <- subset(TK66, subset = percent.mt <50)
TK66 <- subset(TK66, subset = percent.RPS <25)

QUALITY <- TK66@meta.data
write.csv(QUALITY, file = "CSV/E14_Rep2_QC.csv")

VlnPlot(TK66, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK66 <- SCTransform(TK66, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK66 <- RunPCA(TK66, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK66, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK66 <- RunUMAP(TK66, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK66, reduction = "umap", label = F,pt.size = 0.5, 
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

p <- FeaturePlot(TK66, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK66 <- CellCycleScoring(TK66, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK66, group.by = "Phase")

#FeaturePlot(TK66, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK66, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK66 <- FindNeighbors(TK66, dims = 1:20)
TK66 <- FindClusters(TK66, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK66, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK66, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E14_Rep2_DEG-initial.csv")

####1st Cleanup####
#remove 
#12 blood

TK66 <- subset(TK66, idents = c("0","1","2","3","4","5","6","7","8","9",
                                "10","11","13","14","15","16","17","18"))
QUALITY <- TK66@meta.data
write.csv(QUALITY, file = "CSV/E14_Rep2_QC.csv")
TK66 <- SCTransform(TK66, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
TK66 <- RunPCA(TK66, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK66, ndims = 50, reduction = "pca" )
TK66 <- RunUMAP(TK66, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK66, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()
TK66 <- FindNeighbors(TK66, dims = 1:20)
TK66 <- FindClusters(TK66, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK66, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
markers <- FindAllMarkers(TK66, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E14_Rep2_DEG-initial.csv")
#Leave immune cells
#13 blood cells, remove after merging
save(TK66, file = "Robj/E14_Rep2.Robj")


#start from here again
#Merge with HyDD_v1 E14
E14 <- merge(x = TK65, y = list(TK66))
head(E14@meta.data)
E14@meta.data$percent.mt <- NULL
E14@meta.data$percent.RPS <- NULL
E14@meta.data$percent.RPL <- NULL
E14@meta.data$S.Score <- NULL
E14@meta.data$G2M.Score <- NULL
E14@meta.data$Phase <- NULL
E14@meta.data$old.ident <- NULL
E14@meta.data$SCT_snn_res.0.4 <- NULL
E14@meta.data$seurat_clusters <- NULL
save(E14, file = "Robj/E14_Merge.Robj")
Idents(E14) <- "orig.ident"
E14[["percent.sex"]] <- PercentageFeatureSet(E14,
                                             features = c("Xist","Malat1","Tsix"))
E14[["percent.mt"]] <- PercentageFeatureSet(E14, pattern = "^mt-")
E14[["percent.RPS"]] <- PercentageFeatureSet(E14, pattern = "^Rps")
E14[["percent.RPL"]] <- PercentageFeatureSet(E14, pattern = "^Rpl")


VlnPlot(E14, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                          "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)

table(E14@meta.data$orig.ident)
Idents(E14) <- "Age"
E14<- SCTransform(E14, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA

#Run PCA
E14<- RunPCA(E14, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(E14, ndims = 50, reduction = "pca" )
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
E14 <- CellCycleScoring(E14, s.features = m.s.genes$mouseGene,
                        g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)

####harmony####
E14<- RunHarmony(E14, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E14<- RunUMAP(E14, reduction = "harmony", dims = 1:20,
              n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(E14, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(E14, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
DimPlot(E14, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(E14, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
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
           "Avp","Crh","Oxt","Th","Vip","Trh","Rfrp")
p <- FeaturePlot(E14, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

save(E14, file = "Robj/E14_Merge.Robj")



####Cluster-identity####
E14 <- FindNeighbors(E14, dims = 1:20, reduction = "harmony")
E14 <- FindClusters(E14, resolution = 2) #SCT_snn_res.0.8
DimPlot(E14, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
table(E14@active.ident)

#initial markers
markers <- FindAllMarkers(E14, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/E14-Cluster-DEG-initial.csv")


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

p <- FeaturePlot(E14, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = F, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)



#Rename
#LH observed
#ID/Ant ID shown
# POA/SCN maybe?

#0 PMN
#1 SMN
#2 MMN
#3 Thalamus?
#4 PrethaLAmus
#5 SMN
#6 PRETHALAMUS
#7 Tuberal (ARC/VMH)
#8 SMN
#9 Tuberal (PMN)
#10 PVH/SON
#11 MMN
#12 Ant ID?
#13 PRETHALAMUS (HAL EMIENENCE)
#14 NPC
#15 PVH/SON
#16 Foxg1/PVH/SON
#17 Check (Foxg1/MGE?)
#18 neural pro (ascl1)
#19 NPC
#20 MMN
#21 PVH&SON
#22 NPC
#23 neural pro (Neurog2)
#24 SMN
#25 Check (DMH?)
#26 Tuberal (ARC/VMH)
#27 Foxg1 NPC
#28 NPC 
#29 NPC
#30 JUNK
#31 ID
#32 JUNK
#33 NPC
#34 PRETHALAMUS
#35 SMN
#36 ANT ID?
#37 ENDOTHELIAL
#38 IMMUNE
#39 Tuberal (PMN-LH)
#40 Check (PMN?)
#41 Foxg1 NPC
#42 Tuberal (LH)
#43 Immune cell
#44 Immune cell


E14 <- RenameIdents(E14, "0"  = "Tuberal (PMN)", "1" = "SMN", "2" = "MMN", "3" = "Check (Thalamus)",
                    "4" = "Prethalamus", "5" = "SMN", "6" = "Prethalamus", "7" = "Tuberal (ARC/VMH)",
                    "8" = "SMN", "9" = "Tuberal (PMN)", "10" = "PVH/SON", "11" = "MMN", "12" = "Ant ID?",
                    "13" = "Prethalamus", "14" = "NPC", "15" = "PVH/SON", "16" = "Foxg1/PVH/SON",
                    "17" = "Check (MGE?)", "18" = "Neural Pro (Ascl1)", "19" = "NPC", "20"= "MMN",
                    "21" = "PVH/SON", "22" = "NPC", "23" = "Neural Pro (Neurog2)",  "24" = "SMN",
                    "25" = "Check (DMH?)", "26" = "Tuberal (ARC/VMH)", "27" = "NPC (Foxg1)", "28" = "NPC",
                    "29" = "NPC", "30" = "Junk", "31" = "ID?", "32" = "Junk", "33" = "NPC", "34" = "Prethalamus",
                    "35" = "SMN", "36" = "Ant ID?", "37" = "Endothelial", "38" = "Immune cells", 
                    "39" = "Tuberal (PMN-LH)", "40" = "Check (PMN?)", "41" = "NPC (Foxg1)", "42" = "Tuberal (LH)",
                    "43" = "Immune cells", "44" = "Immune cells")

Group <- sort(levels(E14@active.ident))
E14 <- AddMetaData(E14, E14@active.ident, "Cluster_Pass1")
E14@active.ident <- factor(x = E14@active.ident, levels = Group)
E14@meta.data$Cluster_Pass1 <- factor(x = E14@meta.data$Cluster_Pass1, levels = Group)
DimPlot(E14, reduction = "umap", label = F, pt.size = 0.5)

save(E14, file = "Robj/E14_Merge.Robj")
