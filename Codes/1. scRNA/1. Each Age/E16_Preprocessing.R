library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("D:/")
####HyDD2####
#E16.5 
#Fresh
####TK49####
#New v3.1 10x HyDD2
#TK49 - first replicate
TK49<- Read10X_h5("scRNA/TK49/outs/raw_feature_bc_matrix.h5")
#adding column names to avoid barcode clash
colnames(TK49) = paste0("TK49_", colnames(TK49))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK49 <- CreateSeuratObject(counts = TK49, project = "E16",
                           min.cells = 5, min.features = 1000) 
TK49 <- subset(TK49, subset = nCount_RNA > 2000)  

###merge####
TK49 <- RenameIdents(TK49, "TK49"= "E16_Rep1")
TK49 <- AddMetaData(TK49, TK49@active.ident, "Age")
head(TK49@meta.data)
Idents(TK49) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK49[["percent.mt"]] <- PercentageFeatureSet(TK49, pattern = "^mt-")
TK49[["percent.RPS"]] <- PercentageFeatureSet(TK49, pattern = "^Rps")
TK49[["percent.RPL"]] <- PercentageFeatureSet(TK49, pattern = "^Rpl")


#filter some rough values
TK49 <- subset(TK49, subset = percent.mt <50)
TK49 <- subset(TK49, subset = percent.RPS <25)

QUALITY <- TK49@meta.data
write.csv(QUALITY, file = "CSV/E16_Rep1_QC.csv")

VlnPlot(TK49, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK49 <- SCTransform(TK49, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK49 <- RunPCA(TK49, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK49, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK49 <- RunUMAP(TK49, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK49, reduction = "umap", label = F,pt.size = 0.5, 
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
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
m.s.genes$mouseGene
m.g2m.genes$mouseGene

p <- FeaturePlot(TK49, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK49 <- CellCycleScoring(TK49, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK49, group.by = "Phase")

#FeaturePlot(TK49, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK49, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK49 <- FindNeighbors(TK49, dims = 1:20)
TK49 <- FindClusters(TK49, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK49, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK49, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E16_Rep1_DEG-initial.csv")

save(TK49, file = "Robj/E16_Rep1.Robj")

####TK50####
#New v3.1 10x HyDD2
#Second repliate
TK50<- Read10X_h5("scRNA/TK50/outs/raw_feature_bc_matrix.h5") 
#adding column names to avoid barcode clash
colnames(TK50) = paste0("TK50_", colnames(TK50))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK50 <- CreateSeuratObject(counts = TK50, project = "E16",
                           min.cells = 5, min.features = 1000) 
TK50 <- subset(TK50, subset = nCount_RNA > 2000) 

###merge####
TK50 <- RenameIdents(TK50, "TK50"= "E16_Rep2")
TK50 <- AddMetaData(TK50, TK50@active.ident, "Age")
head(TK50@meta.data)
Idents(TK50) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK50[["percent.mt"]] <- PercentageFeatureSet(TK50, pattern = "^mt-")
TK50[["percent.RPS"]] <- PercentageFeatureSet(TK50, pattern = "^Rps")
TK50[["percent.RPL"]] <- PercentageFeatureSet(TK50, pattern = "^Rpl")

#filter some rough values
TK50 <- subset(TK50, subset = percent.mt <50)
TK50 <- subset(TK50, subset = percent.RPS <25)

QUALITY <- TK50@meta.data
write.csv(QUALITY, file = "CSV/E16_Rep2_QC.csv")

VlnPlot(TK50, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK50 <- SCTransform(TK50, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK50 <- RunPCA(TK50, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK50, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK50 <- RunUMAP(TK50, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK50, reduction = "umap", label = F,pt.size = 0.5, 
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

p <- FeaturePlot(TK50, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK50 <- CellCycleScoring(TK50, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK50, group.by = "Phase")

#FeaturePlot(TK50, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK50, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK50 <- FindNeighbors(TK50, dims = 1:20)
TK50 <- FindClusters(TK50, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK50, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK50, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E16_Rep2_DEG-initial.csv")

#Leave immune cells
save(TK50, file = "Robj/E16_Rep2.Robj")


#start from here again
#Merge with HyDD_v1 E16
E16 <- merge(x = TK49, y = list(TK50))
head(E16@meta.data)
E16@meta.data$percent.mt <- NULL
E16@meta.data$percent.RPS <- NULL
E16@meta.data$percent.RPL <- NULL
E16@meta.data$S.Score <- NULL
E16@meta.data$G2M.Score <- NULL
E16@meta.data$Phase <- NULL
E16@meta.data$old.ident <- NULL
E16@meta.data$SCT_snn_res.0.4 <- NULL
E16@meta.data$seurat_clusters <- NULL
save(E16, file = "Robj/E16_Merge.Robj")
load(file = "Robj/E16_Merge.Robj")
Idents(E16) <- "orig.ident"
E16[["percent.sex"]] <- PercentageFeatureSet(E16,
                                             features = c("Xist","Malat1","Tsix"))
E16[["percent.mt"]] <- PercentageFeatureSet(E16, pattern = "^mt-")
E16[["percent.RPS"]] <- PercentageFeatureSet(E16, pattern = "^Rps")
E16[["percent.RPL"]] <- PercentageFeatureSet(E16, pattern = "^Rpl")


VlnPlot(E16, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                          "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)

table(E16@meta.data$orig.ident)
Idents(E16) <- "Age"
E16<- SCTransform(E16, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA

#Run PCA
E16<- RunPCA(E16, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(E16, ndims = 50, reduction = "pca" )
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
E16 <- CellCycleScoring(E16, s.features = m.s.genes$mouseGene,
                        g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)

####harmony####
E16<- RunHarmony(E16, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E16<- RunUMAP(E16, reduction = "harmony", dims = 1:20,
              n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(E16, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(E16, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
DimPlot(E16, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(E16, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
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
p <- FeaturePlot(E16, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

save(E16, file = "Robj/E16_Merge.Robj")

####Cluster-identity####
E16 <- FindNeighbors(E16, dims = 1:20, reduction = "harmony")
E16 <- FindClusters(E16, resolution = 2) #SCT_snn_res.0.8
DimPlot(E16, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(E16, file = "Robj/E16_Merge.Robj")

#initial markers
markers <- FindAllMarkers(E16, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/E16-Cluster-DEG-initial.csv")





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

p <- FeaturePlot(E16, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = F, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)



#Rename
#E16
#also observed similar in E15 but maybe rise of regional progenitors
#PMN shift?
#0 SMN Tac1
#1 ARC/VMH
#2 Junk?
#3 PMN
#4 Ant ID?
#5 Pomc/Npy/Sst (Arc?)
#6 Bsx/Cck (PMN?)
#7 Cck (MMN)
#8 PMN?
#9 Glial progenitor - Tanycyte?
#10 Cck (MMN)
#11 Glial progenitor
#12 Glial progenitor
#13 Nts/Tac1/Cck (MMN)
#14 Prethal/ID
#15 Glial progenitor
#16 MMN
#17 SMN? is pvn/son there
#18 Trh/Pnoc (POA/SCN?)
#19 MMN
#20 Penk (POA/SCN?)
#21 Gal (PMN)
#22 Junk
#23 Immune cells? Glial progenitor?
#24 SMN?
#25 ARC/VMH
#26 Prethal/ID (glial?)
#27 Rora/Rorb (POA/SCN?)
#28 SMN?  is pvn/son there
#29 Glial progenitor? - OPC/Oligo?
#30 Glial progenitor?
#31 Glial progenitor?
#32 Immune cells?
#33 Avp/Gal/Oxt (PVN?)
#34 Glial progenitor?
#35 Junk
#36 Glial progenitor?
#37 Immune cells?
#38 Immune cells
#39 Immune cells
#40 Hcrt/Npvf
#41 Glial progenitor?

E16 <- RenameIdents(E16, "0" = "Tac1 (SMN)","1" = "ARC/VMH", "2" = "Junk", "3" = "PMN",
                    "4" = "Check (Ant ID?)", "5" = "Pomc/Npy/Sst (Arc)", "6" = "Bsx/Cck (PMN?)",
                    "7" = "Cck (MMN)", "8" = "Check (PMN?)", "9" = "Glial Pro? (Tanycyte?)",
                    "10" = "Cck (MMN)", "11" = "Glial Pro?", "12" = "Glial Pro?", "13"  = "Nts/Tac1/Cck (MMN)",
                    "14" = "PreThal/ID", "15" = "Glial Pro?", "16" = "MMN", "17" = "Check (SMN/PVN/SON?)",
                    "18" = "Trh/Pnoc (POA/SCN?)", "19" = "MMN", "20" = "Penk (POA/SCN?)", "21" = "Gal (PMN)",
                    "22"= "Junk", "23" = "Glial Pro?", "24" = "Check (SMN?)", "25" = "ARC/VMH",
                    "26" = "PreThal/ID (Glial?)", "27" = "Rora/Rorb (POA/SCN?)", "28" = "Check (SMN/PVN/SON?)",
                    "29"  = "Glial Pro?", "30"  = "Glial Pro?", "31"  = "Glial Pro?", "32" = "Immune cells?",
                    "33" = "Avp/Gal/Oxt (PVN/SON?)", "34" = "Glial Pro?", "35" = "Junk", "36" = "Glial Pro?",
                    "37" = "Immune cells?", "38" = "Immune cells?", "39" = "Immune cells?", "40" = "Hcrt/Npvf",
                    "41" = "Glial Pro?"  )
Group <- sort(levels(E16@active.ident))
E16 <- AddMetaData(E16, E16@active.ident, "Cluster_Pass1")
E16@active.ident <- factor(x = E16@active.ident, levels = Group)
E16@meta.data$Cluster_Pass1 <- factor(x = E16@meta.data$Cluster_Pass1, levels = Group)
DimPlot(E16, reduction = "umap", label = F, pt.size = 0.5)

save(E16, file = "Robj/E16_Merge.Robj")

