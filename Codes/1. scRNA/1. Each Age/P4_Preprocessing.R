library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("D:/")
####HyDD2####
#P4.5 
#Fresh
####TK82####
#New v3.1 10x HyDD2
#TK82 - first replicate
TK82<- Read10X_h5("scRNA/TK82/outs/raw_feature_bc_matrix.h5") #500 million reads/sample
#adding column names to avoid barcode clash
colnames(TK82) = paste0("TK82_", colnames(TK82))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK82 <- CreateSeuratObject(counts = TK82, project = "P4",
                           min.cells = 5, min.features = 1000) #1000 if depth is higher

#Filter based on these units.  #not by mitochondria and ribosomal -> aging sample
TK82 <- subset(TK82, subset = nCount_RNA > 2000)  #UMI #1000/2000/3000

###merge####
TK82 <- RenameIdents(TK82, "TK82"= "P4_Rep1")
TK82 <- AddMetaData(TK82, TK82@active.ident, "Age")
head(TK82@meta.data)
Idents(TK82) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK82[["percent.mt"]] <- PercentageFeatureSet(TK82, pattern = "^mt-")
TK82[["percent.RPS"]] <- PercentageFeatureSet(TK82, pattern = "^Rps")
TK82[["percent.RPL"]] <- PercentageFeatureSet(TK82, pattern = "^Rpl")
#mito.genes <- grep(pattern = "^mt-", x = rownames(x = TK82@assays$SCT@data), value = TRUE)
#percent.mito <- colSums(expm1(TK82@assays$SCT@data[mito.genes, ]))/colSums(expm1(TK82@assays$SCT@data))
#TK82 <- AddMetaData(TK82, percent.mito, "percent.mito")

#filter some rough values
TK82 <- subset(TK82, subset = percent.mt <50)
TK82 <- subset(TK82, subset = percent.RPS <25)

QUALITY <- TK82@meta.data
write.csv(QUALITY, file = "CSV/P4_Rep1_QC.csv")

VlnPlot(TK82, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK82 <- SCTransform(TK82, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK82 <- RunPCA(TK82, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK82, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK82 <- RunUMAP(TK82, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK82, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()

####Check####
genes = c("Foxd1","Foxg1","Shh","Fst","Pax6","SP4","Sp9","Nkx2-1","Rax","Foxa1","Foxa2")
genes <- c("Dlx2", "Olig2", "Nkx2-2", "Olfml3", "Chrd", "Nov", "Car2",
           "Sst", "Fgf18", "Six6", "Dbx1", "Bmp2", "Bmp7", "Vax1", "Vax2", "Emx2", "Nkx2-4", "Sfrp2", "Pitx2")
genes <- c("Six3", "Lef1", "Axin2", "Top2a", "Cdc20", "Ccnb1", "Ccne1")
genes <- c("Sim1", "Otp", "Emx2", "Lhx5", "Lhx1", "Isl1", "Nr5a1", "Wnt8a", "Wnt8b", "Fgf10", "Pax2", "Vax1","Vax2")
genes <- c("Xist","Malat1","Tsix","Rps4x","Ddx3x","Ddx3y")
genes <- c("Avp", "Spp1", "Cst3", "Rorb")
genes <- c("Lhx2", "Lhx9", "Pitx3", "Gbx2", "Tcf7l2", "Gata2", "Gata3", "Hmx2","Hmx3","Nr5a1")
genes <- c("Foxd1","Lhx8","Shh","Fst","Pax6","Sp8","Sp9","Nkx2-1","Rax","Foxa1","Foxa2","Lhx9","Sim1","Otp","Nr5a1",
           "Agrp","Npy","Gal","Slc32a1","Slc17a6","Pomc","Cartpt","Hcrt","Nts","Sst","Ghrh","Gnrh1","Lepr",
           "Avp","Crh","Oxt","Th","Vip","Trh","Npvf")
genes <- c("Aldh1a1","Ntsr2","Gfap","Slc1a3","Rax","Foxj1","Mobp","Mog","Olig1","Ptgds","Cd74")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
m.s.genes$mouseGene
m.g2m.genes$mouseGene

p <- FeaturePlot(TK82, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK82 <- CellCycleScoring(TK82, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK82, group.by = "Phase")

#FeaturePlot(TK82, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK82, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK82 <- FindNeighbors(TK82, dims = 1:20)
TK82 <- FindClusters(TK82, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK82, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK82, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/P4_Rep1_DEG-initial.csv")
save(TK82, file = "Robj/P4_Rep1.Robj")

####TK83####
#New v3.1 10x HyDD2
#Second repliate
TK83<- Read10X_h5("scRNA/TK83/outs/raw_feature_bc_matrix.h5") #500 million reads/sample
#adding column names to avoid barcode clash
colnames(TK83) = paste0("TK83_", colnames(TK83))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK83 <- CreateSeuratObject(counts = TK83, project = "P4",
                           min.cells = 5, min.features = 1000) #1000 if depth is higher

#Filter based on these units.  #not by mitochondria and ribosomal -> aging sample
TK83 <- subset(TK83, subset = nCount_RNA > 2000)  #UMI #1000/2000/3000

###merge####
TK83 <- RenameIdents(TK83, "TK83"= "P4_Rep2")
TK83 <- AddMetaData(TK83, TK83@active.ident, "Age")
head(TK83@meta.data)
Idents(TK83) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK83[["percent.mt"]] <- PercentageFeatureSet(TK83, pattern = "^mt-")
TK83[["percent.RPS"]] <- PercentageFeatureSet(TK83, pattern = "^Rps")
TK83[["percent.RPL"]] <- PercentageFeatureSet(TK83, pattern = "^Rpl")
#mito.genes <- grep(pattern = "^mt-", x = rownames(x = TK83@assays$SCT@data), value = TRUE)
#percent.mito <- colSums(expm1(TK83@assays$SCT@data[mito.genes, ]))/colSums(expm1(TK83@assays$SCT@data))
#TK83 <- AddMetaData(TK83, percent.mito, "percent.mito")

#filter some rough values
TK83 <- subset(TK83, subset = percent.mt <50)
TK83 <- subset(TK83, subset = percent.RPS <25)

QUALITY <- TK83@meta.data
write.csv(QUALITY, file = "CSV/P4_Rep2_QC.csv")

VlnPlot(TK83, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK83 <- SCTransform(TK83, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK83 <- RunPCA(TK83, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK83, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK83 <- RunUMAP(TK83, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK83, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()

####Check####
genes = c("Foxd1","Foxg1","Shh","Fst","Pax6","SP4","Sp9","Nkx2-1","Rax","Foxa1","Foxa2")
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

p <- FeaturePlot(TK83, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK83 <- CellCycleScoring(TK83, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK83, group.by = "Phase")

#FeaturePlot(TK83, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK83, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK83 <- FindNeighbors(TK83, dims = 1:20)
TK83 <- FindClusters(TK83, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK83, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK83, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/P4_Rep2_DEG-initial.csv")

#Leave immune cells
save(TK83, file = "Robj/P4_Rep2.Robj")

####TK84####
#New v3.1 10x HyDD2
#TK84 - first replicate
TK84<- Read10X_h5("scRNA/TK84/outs/raw_feature_bc_matrix.h5") #500 million reads/sample
#adding column names to avoid barcode clash
colnames(TK84) = paste0("TK84_", colnames(TK84))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK84 <- CreateSeuratObject(counts = TK84, project = "P4",
                           min.cells = 5, min.features = 1000) #1000 if depth is higher

#Filter based on these units.  #not by mitochondria and ribosomal -> aging sample
TK84 <- subset(TK84, subset = nCount_RNA > 2000)  #UMI #1000/2000/3000

###merge####
TK84 <- RenameIdents(TK84, "TK84"= "P4_Rep3")
TK84 <- AddMetaData(TK84, TK84@active.ident, "Age")
head(TK84@meta.data)
Idents(TK84) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK84[["percent.mt"]] <- PercentageFeatureSet(TK84, pattern = "^mt-")
TK84[["percent.RPS"]] <- PercentageFeatureSet(TK84, pattern = "^Rps")
TK84[["percent.RPL"]] <- PercentageFeatureSet(TK84, pattern = "^Rpl")
#mito.genes <- grep(pattern = "^mt-", x = rownames(x = TK84@assays$SCT@data), value = TRUE)
#percent.mito <- colSums(expm1(TK84@assays$SCT@data[mito.genes, ]))/colSums(expm1(TK84@assays$SCT@data))
#TK84 <- AddMetaData(TK84, percent.mito, "percent.mito")

#filter some rough values
TK84 <- subset(TK84, subset = percent.mt <50)
TK84 <- subset(TK84, subset = percent.RPS <25)

QUALITY <- TK84@meta.data
write.csv(QUALITY, file = "CSV/P4_Rep3_QC.csv")

VlnPlot(TK84, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK84 <- SCTransform(TK84, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK84 <- RunPCA(TK84, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK84, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK84 <- RunUMAP(TK84, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK84, reduction = "umap", label = F,pt.size = 0.5, 
        cols = , ncol = 2) + NoLegend()

####Check####
genes = c("Foxd1","Foxg1","Shh","Fst","Pax6","SP4","Sp9","Nkx2-1","Rax","Foxa1","Foxa2")
genes <- c("Dlx2", "Olig2", "Nkx2-2", "Olfml3", "Chrd", "Nov", "Car2",
           "Sst", "Fgf18", "Six6", "Dbx1", "Bmp2", "Bmp7", "Vax1", "Vax2", "Emx2", "Nkx2-4", "Sfrp2", "Pitx2")
genes <- c("Six3", "Lef1", "Axin2", "Top2a", "Cdc20", "Ccnb1", "Ccne1")
genes <- c("Sim1", "Otp", "Emx2", "Lhx5", "Lhx1", "Isl1", "Nr5a1", "Wnt8a", "Wnt8b", "Fgf10", "Pax2", "Vax1","Vax2")
genes <- c("Xist","Malat1","Tsix","Rps4x","Ddx3x","Ddx3y")
genes <- c("Avp", "Spp1", "Cst3", "Rorb")
genes <- c("Lhx2", "Lhx9", "Pitx3", "Gbx2", "Tcf7l2", "Gata2", "Gata3", "Hmx2","Hmx3","Nr5a1")
genes <- c("Foxd1","Lhx8","Shh","Fst","Pax6","SP4","Sp9","Nkx2-1","Rax","Foxa1","Foxa2","Lhx9","Sim1","Otp","Nr5a1",
           "Agrp","Npy","Gal","Slc32a1","Slc17a6","Pomc","Cartpt","Hcrt","Nts","Sst","Ghrh","Gnrh1","Lepr",
           "Avp","Crh","Oxt","Th","Vip","Trh","Npvf")
genes <- c("Aldh1a1","Ntsr2","Gfap","Slc1a3","Rax","Foxj1","Mobp","Mog","Olig1","Ptgds","Cd74")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
m.s.genes$mouseGene
m.g2m.genes$mouseGene

p <- FeaturePlot(TK84, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK84 <- CellCycleScoring(TK84, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK84, group.by = "Phase")

#FeaturePlot(TK84, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK84, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK84 <- FindNeighbors(TK84, dims = 1:20)
TK84 <- FindClusters(TK84, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK84, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK84, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/P4_Rep3_DEG-initial.csv")
save(TK84, file = "Robj/P4_Rep3.Robj")

#start from here again
#Merge with HyDD_v1 P4
P4 <- merge(x = TK82, y = list(TK83,TK84))
head(P4@meta.data)
P4@meta.data$percent.mt <- NULL
P4@meta.data$percent.RPS <- NULL
P4@meta.data$percent.RPL <- NULL
P4@meta.data$S.Score <- NULL
P4@meta.data$G2M.Score <- NULL
P4@meta.data$Phase <- NULL
P4@meta.data$old.ident <- NULL
P4@meta.data$SCT_snn_res.0.4 <- NULL
P4@meta.data$seurat_clusters <- NULL
save(P4, file = "Robj/P4_Merge.Robj")
load(file = "Robj/P4_Merge.Robj")
Idents(P4) <- "orig.ident"
P4[["percent.sex"]] <- PercentageFeatureSet(P4,
                                            features = c("Xist","Malat1","Tsix"))
P4[["percent.mt"]] <- PercentageFeatureSet(P4, pattern = "^mt-")
P4[["percent.RPS"]] <- PercentageFeatureSet(P4, pattern = "^Rps")
P4[["percent.RPL"]] <- PercentageFeatureSet(P4, pattern = "^Rpl")


VlnPlot(P4, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                         "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)

table(P4@meta.data$orig.ident)
Idents(P4) <- "Age"
P4<- SCTransform(P4, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA

#Run PCA
P4<- RunPCA(P4, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(P4, ndims = 50, reduction = "pca" )
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
P4 <- CellCycleScoring(P4, s.features = m.s.genes$mouseGene,
                       g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)

####harmony####
P4<- RunHarmony(P4, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
P4<- RunUMAP(P4, reduction = "harmony", dims = 1:20,
             n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(P4, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(P4, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
DimPlot(P4, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(P4, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
        split.by = "Phase", ncol = 2) + NoLegend()

####Check####
genes = c("Foxd1","Foxg1","Shh","Fst","Pax6","SP4","Sp9","Nkx2-1","Rax","Foxa1","Foxa2")
genes <- c("Dlx2", "Olig2", "Nkx2-2", "Olfml3", "Chrd", "Nov", "Car2",
           "Sst", "Fgf18", "Six6", "Dbx1", "Bmp2", "Bmp7", "Vax1", "Vax2", "Emx2", "Nkx2-4", "Sfrp2", "Pitx2")
genes <- c("Six3", "Lef1", "Axin2", "Top2a", "Cdc20", "Ccnb1", "Ccne1")
genes <- c("Sim1", "Otp", "Emx2", "Lhx5", "Lhx1", "Isl1", "Nr5a1", "Wnt8a", "Wnt8b", "Fgf10", "Pax2", "Vax1","Vax2")
genes <- c("Xist","Malat1","Tsix","Rps4x","Ddx3x","Ddx3y")
genes <- c("Avp", "Spp1", "Cst3", "Rorb")
genes <- c("Lhx2", "Lhx9", "Pitx3", "Gbx2", "Tcf7l2", "Gata2", "Gata3", "Hmx2","Hmx3")
genes <- c("Foxd1","Lhx8","Shh","Fst","Pax6","Sp8","Sp9","Nkx2-1","Rax","Foxa1","Foxa2","Lhx9","Sim1","Otp",
           "Agrp","Npy","Gal","Slc32a1","Slc17a6","Pomc","Cartpt","Hcrt","Nts","Sst","Ghrh","Gnrh1","Lepr",
           "Avp","Crh","Oxt","Th","Vip","Trh","Npvf")
p <- FeaturePlot(P4, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

save(P4, file = "Robj/P4_Merge.Robj")

####Cluster-identity####
P4 <- FindNeighbors(P4, dims = 1:20, reduction = "harmony")
P4 <- FindClusters(P4, resolution = 2) #SCT_snn_res.0.8
DimPlot(P4, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(P4, file = "Robj/P4_Merge.Robj")

#initial markers
markers <- FindAllMarkers(P4, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/P4-Cluster-DEG-initial.csv")





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

p <- FeaturePlot(P4, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = F, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)



#Rename
#P4
P4
#0 neuron
#1 glia
#2 glia 
#3 glia
#4 glia
#5 neuron
#6 neuron
#7 glia
#8 glia
#9 neuron
#10 immune
#11 glia
#12 neuron
#13 glia
#14 immune
#15 neuron (remove?)
#16 neuron
#17 neuron
#18 glia
#19 glia
#20 glia
#21 glia
#22 neuron
#23 glia
#24 immune
#25 junk?
#26 glia
#27 neuron
#28 neuron
#29 neuron
#30 endothelial
#31 junk
#32 glia
#33 glia
#34 endothelial
#35 pericyte?
#36 bam
#37 immune
#38 glia
#39 endothelial
#40 junk
#41 msc
#42 neuron
#43 neuron
#44 neuron
#45 blood
#46 glia
#47 immune
#48 pituitary
#49 pituitary
#50 glia
#51 blood
#52 neuron

names <- levels(P4@active.ident)
P4 <- RenameIdents(P4, "0" = "Neuron", "1"= "Glia", "2"= "Glia", "3"= "Glia", "4"= "Glia",
                   "5"= "Neuron", "6"= "Neuron",
                   "7"= "Glia", "8"= "Glia", "9"= "Neuron", "10"= "Immune", "11"= "Glia", 
                   "12"= "Neuron", "13"= "Glia", "14"= "Immune", "15" = "Neuron", "16"= "Neuron", "17"= "Neuron",
                   "18"= "Glia", "19"= "Glia", "20"= "Glia", "21"= "Glia", "22"= "Neuron", 
                   "23"= "Glia", "24"= "Immune", "25"= "Junk", "26"= "Glia",
                   "27"= "Neuron", "28"= "Neuron", "29"= "Neuron",
                   "30"="Endothelial", "31"= "Junk", "32"= "Glia", "33"= "Glia", 
                   "34"="Endothelial", "35"= "Immune", "36"= "Immune", "37"= "Immune",
                   "38"= "Glia", "39"="Endothelial", "40"= "Junk", "41"= "MSC",
                   "42"= "Neuron", "43"= "Neuron", "44"= "Neuron", 
                   "45"="Blood", "46"= "Glia", "47"= "Immune", "48"= "Pituitary", "49"= "Pituitary",
                   "50"= "Glia", "51"="Blood", "52" = "Neuron"  )

DimPlot(P4, reduction = "umap", label = F,pt.size = 0.1)
Group <- sort(levels(P4@active.ident))
P4 <- AddMetaData(P4, P4@active.ident, "Cluster_Pass1")
P4@active.ident <- factor(x = P4@active.ident, levels = Group)
P4@meta.data$Cluster_Pass1 <- factor(x = P4@meta.data$Cluster_Pass1, levels = Group)
DimPlot(P4, reduction = "umap", label = F, pt.size = 0.5)
save(P4, file = "Robj/P4_Merge.Robj")

