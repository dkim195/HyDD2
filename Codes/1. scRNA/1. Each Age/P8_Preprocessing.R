library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("D:/")
####HyDD2####
#P8.5 
#Fresh
####TK85####
#New v3.1 10x HyDD2
#TK85 - first replicate
TK85<- Read10X_h5("scRNA/TK85/outs/raw_feature_bc_matrix.h5") #500 million reads/sample
#adding column names to avoid barcode clash
colnames(TK85) = paste0("TK85_", colnames(TK85))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK85 <- CreateSeuratObject(counts = TK85, project = "P8",
                           min.cells = 5, min.features = 1000) #1000 if depth is higher

#Filter based on these units.  #not by mitochondria and ribosomal -> aging sample
TK85 <- subset(TK85, subset = nCount_RNA > 2000)  #UMI #1000/2000/3000

###merge####
TK85 <- RenameIdents(TK85, "TK85"= "P8_Rep1")
TK85 <- AddMetaData(TK85, TK85@active.ident, "Age")
head(TK85@meta.data)
Idents(TK85) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK85[["percent.mt"]] <- PercentageFeatureSet(TK85, pattern = "^mt-")
TK85[["percent.RPS"]] <- PercentageFeatureSet(TK85, pattern = "^Rps")
TK85[["percent.RPL"]] <- PercentageFeatureSet(TK85, pattern = "^Rpl")
#mito.genes <- grep(pattern = "^mt-", x = rownames(x = TK85@assays$SCT@data), value = TRUE)
#percent.mito <- colSums(expm1(TK85@assays$SCT@data[mito.genes, ]))/colSums(expm1(TK85@assays$SCT@data))
#TK85 <- AddMetaData(TK85, percent.mito, "percent.mito")

#filter some rough values
TK85 <- subset(TK85, subset = percent.mt <50)
TK85 <- subset(TK85, subset = percent.RPS <25)

QUALITY <- TK85@meta.data
write.csv(QUALITY, file = "CSV/P8_Rep1_QC.csv")

VlnPlot(TK85, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK85 <- SCTransform(TK85, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK85 <- RunPCA(TK85, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK85, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK85 <- RunUMAP(TK85, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK85, reduction = "umap", label = F,pt.size = 0.5, 
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

p <- FeaturePlot(TK85, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK85 <- CellCycleScoring(TK85, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK85, group.by = "Phase")

#FeaturePlot(TK85, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK85, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK85 <- FindNeighbors(TK85, dims = 1:20)
TK85 <- FindClusters(TK85, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK85, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK85, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/P8_Rep1_DEG-initial.csv")
save(TK85, file = "Robj/P8_Rep1.Robj")

####TK90####
#New v3.1 10x HyDD2
#Second repliate
TK90<- Read10X_h5("scRNA/TK90/outs/raw_feature_bc_matrix.h5") #500 million reads/sample
#adding column names to avoid barcode clash
colnames(TK90) = paste0("TK90_", colnames(TK90))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK90 <- CreateSeuratObject(counts = TK90, project = "P8",
                           min.cells = 5, min.features = 1000) #1000 if depth is higher

#Filter based on these units.  #not by mitochondria and ribosomal -> aging sample
TK90 <- subset(TK90, subset = nCount_RNA > 2000)  #UMI #1000/2000/3000

###merge####
TK90 <- RenameIdents(TK90, "TK90"= "P8_Rep2")
TK90 <- AddMetaData(TK90, TK90@active.ident, "Age")
head(TK90@meta.data)
Idents(TK90) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK90[["percent.mt"]] <- PercentageFeatureSet(TK90, pattern = "^mt-")
TK90[["percent.RPS"]] <- PercentageFeatureSet(TK90, pattern = "^Rps")
TK90[["percent.RPL"]] <- PercentageFeatureSet(TK90, pattern = "^Rpl")
#mito.genes <- grep(pattern = "^mt-", x = rownames(x = TK90@assays$SCT@data), value = TRUE)
#percent.mito <- colSums(expm1(TK90@assays$SCT@data[mito.genes, ]))/colSums(expm1(TK90@assays$SCT@data))
#TK90 <- AddMetaData(TK90, percent.mito, "percent.mito")

#filter some rough values
TK90 <- subset(TK90, subset = percent.mt <50)
TK90 <- subset(TK90, subset = percent.RPS <25)

QUALITY <- TK90@meta.data
write.csv(QUALITY, file = "CSV/P8_Rep2_QC.csv")

VlnPlot(TK90, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK90 <- SCTransform(TK90, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK90 <- RunPCA(TK90, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK90, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK90 <- RunUMAP(TK90, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK90, reduction = "umap", label = F,pt.size = 0.5, 
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

p <- FeaturePlot(TK90, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK90 <- CellCycleScoring(TK90, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK90, group.by = "Phase")

#FeaturePlot(TK90, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK90, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK90 <- FindNeighbors(TK90, dims = 1:20)
TK90 <- FindClusters(TK90, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK90, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK90, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/P8_Rep2_DEG-initial.csv")

#Leave immune cells
save(TK90, file = "Robj/P8_Rep2.Robj")


#start from here again
#Merge with HyDD_v1 P8
P8 <- merge(x = TK85, y = list(TK90))
head(P8@meta.data)
P8@meta.data$percent.mt <- NULL
P8@meta.data$percent.RPS <- NULL
P8@meta.data$percent.RPL <- NULL
P8@meta.data$S.Score <- NULL
P8@meta.data$G2M.Score <- NULL
P8@meta.data$Phase <- NULL
P8@meta.data$old.ident <- NULL
P8@meta.data$SCT_snn_res.0.4 <- NULL
P8@meta.data$seurat_clusters <- NULL
save(P8, file = "Robj/P8_Merge.Robj")
load(file = "Robj/P8_Merge.Robj")
Idents(P8) <- "orig.ident"
P8[["percent.sex"]] <- PercentageFeatureSet(P8,
                                             features = c("Xist","Malat1","Tsix"))
P8[["percent.mt"]] <- PercentageFeatureSet(P8, pattern = "^mt-")
P8[["percent.RPS"]] <- PercentageFeatureSet(P8, pattern = "^Rps")
P8[["percent.RPL"]] <- PercentageFeatureSet(P8, pattern = "^Rpl")


VlnPlot(P8, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                          "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)

table(P8@meta.data$orig.ident)
Idents(P8) <- "Age"
P8<- SCTransform(P8, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA

#Run PCA
P8<- RunPCA(P8, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(P8, ndims = 50, reduction = "pca" )
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
P8 <- CellCycleScoring(P8, s.features = m.s.genes$mouseGene,
                        g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)

####harmony####
P8<- RunHarmony(P8, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
P8<- RunUMAP(P8, reduction = "harmony", dims = 1:20,
              n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(P8, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(P8, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
DimPlot(P8, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(P8, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
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
           "Avp","Crh","Oxt","Th","Vip","Trh","Npvf")
p <- FeaturePlot(P8, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

save(P8, file = "Robj/P8_Merge.Robj")

####Cluster-identity####
P8 <- FindNeighbors(P8, dims = 1:20, reduction = "harmony")
P8 <- FindClusters(P8, resolution = 2) #SCT_snn_res.0.8
DimPlot(P8, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(P8, file = "Robj/P8_Merge.Robj")

#initial markers
markers <- FindAllMarkers(P8, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/P8-Cluster-DEG-initial.csv")





#genes to plot
genes <- c("Cdkn1c","Gadd45g","Pkib","Dlk1","Hes6","Ass1") #AH Neural Progenitor
genes <- c("D930028M14Rik", "Onecut2", "Sncg", "Snhg11", "Ina", "Cartpt") #Ant ID
genes <- c("Pomc", "Chchd10", "S100a10", "Six6","Rbp1","Sox14") #ARC
genes <- c("Sst","Isl1","Six6","Otp","Dlk1","Cited1") #DMH
genes <- c("Lhx8", "Snhg11","Cadm1","Pcsk1n","Mapt","Gabra2") #ID/TT
genes <- c("Apoe","Lgals1","Igfbp7","Sepp1","Serpinh1","Tyrobp") #Immune cell
genes <- c("Pmch","Lhx9","Chd3os") #LH?
genes <- c("Dlx2","Lhx8","Dlx6os1","Sp9","Arx","Foxg1") #MGE
genes <- c("Lhx1os","PcP8","Foxb1","Lhx1","Nhlh2","Lhx5") #MMN
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

p <- FeaturePlot(P8, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = F, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)



#Rename
P8
#0 neuron
#1 glia
#2 glia
#3 glia
#4 glia
#5 glia - opc cluster?
#6 neuron
#7 neuron
#8 glia
#9 immune
#10 glia astrocyte cluster
#11 glia - epen/tany cluster
#12 endothelial
#13 neuron
#14 glia
#15 immune
#16 glia
#17 glia
#18 neuron
#19 immune
#20 glia
#21 immune
#22 endothelial
#23 neuron
#24 glia
#25 glia
#26 neuron
#27 bam
#28 neuron
#29 neuron
#30 blood
#31 neuron
#32 immune
#33 glia
#34 msc
#35 endothelial
#36 glia
#37 junk
#38 OPC?/junk
#39 neuron
#40 neuron
#41 glia
#42 neuron
#43 pituitary
#44 blood
#45 glia
#46 glia
#47 glia

names <- levels(P8@active.ident)
P8 <- RenameIdents(P8,"0"="Neuron", "1"="Glia", "2"="Glia", "3"="Glia", "4"="Glia", "5"="Glia",
                   "6"="Neuron", "7"="Neuron", "8"="Glia", "9"="Immune", "10"="Glia", "11"="Glia", 
                   "12"="Endothelial", "13"="Neuron", "14"="Glia", "15"="Immune", "16"="Glia", "17"="Glia",
                   "18"="Neuron", "19"="Immune", "20"="Glia", "21"="Immune", "22"="Endothelial", 
                   "23"="Neuron", "24"="Glia", "25"="Glia", "26"="Neuron", "27"="Immune",
                   "28"="Neuron", "29"="Neuron", "30"="Blood",
                   "31"="Neuron", "32"="Immune", "33"="Glia", 
                   "34"="MSC", "35"="Endothelial", "36"="Glia",
                   "37"="Junk", "38"="Junk", "39"="Neuron", "40"="Neuron",
                   "41"="Glia", "42"="Neuron", "43"="Pituitary", "44"="Blood", 
                   "45"="Glia", "46"="Glia", "47"="Glia")


DimPlot(P8, reduction = "umap", label = F,pt.size = 0.1)


Group <- sort(levels(P8@active.ident))
P8 <- AddMetaData(P8, P8@active.ident, "Cluster_Pass1")
P8@active.ident <- factor(x = P8@active.ident, levels = Group)
P8@meta.data$Cluster_Pass1 <- factor(x = P8@meta.data$Cluster_Pass1, levels = Group)
DimPlot(P8, reduction = "umap", label = F, pt.size = 0.5)

save(P8, file = "Robj/P8_Merge.Robj")


