library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("D:/")
####HyDD2####
#E18.5 
#Fresh
####TK69####
#New v3.1 10x HyDD2
#TK69 - first replicate
TK69<- Read10X_h5("scRNA/TK69/outs/raw_feature_bc_matrix.h5") 
#adding column names to avoid barcode clash
colnames(TK69) = paste0("TK69_", colnames(TK69))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK69 <- CreateSeuratObject(counts = TK69, project = "E18",
                           min.cells = 5, min.features = 1000)

TK69 <- subset(TK69, subset = nCount_RNA > 2000)

###merge####
TK69 <- RenameIdents(TK69, "TK69"= "E18_Rep1")
TK69 <- AddMetaData(TK69, TK69@active.ident, "Age")
head(TK69@meta.data)
Idents(TK69) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK69[["percent.mt"]] <- PercentageFeatureSet(TK69, pattern = "^mt-")
TK69[["percent.RPS"]] <- PercentageFeatureSet(TK69, pattern = "^Rps")
TK69[["percent.RPL"]] <- PercentageFeatureSet(TK69, pattern = "^Rpl")

#filter some rough values
TK69 <- subset(TK69, subset = percent.mt <50)
TK69 <- subset(TK69, subset = percent.RPS <25)

QUALITY <- TK69@meta.data
write.csv(QUALITY, file = "CSV/E18_Rep1_QC.csv")

VlnPlot(TK69, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK69 <- SCTransform(TK69, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK69 <- RunPCA(TK69, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK69, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK69 <- RunUMAP(TK69, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK69, reduction = "umap", label = F,pt.size = 0.5, 
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

p <- FeaturePlot(TK69, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK69 <- CellCycleScoring(TK69, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK69, group.by = "Phase")

#FeaturePlot(TK69, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK69, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK69 <- FindNeighbors(TK69, dims = 1:20)
TK69 <- FindClusters(TK69, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK69, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK69, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E18_Rep1_DEG-initial.csv")
save(TK69, file = "Robj/E18_Rep1.Robj")

####TK70####
#New v3.1 10x HyDD2
#Second repliate
TK70<- Read10X_h5("scRNA/TK70/outs/raw_feature_bc_matrix.h5")
#adding column names to avoid barcode clash
colnames(TK70) = paste0("TK70_", colnames(TK70))
#Create Seurat object. Cutoff 200 genes in 5 cells
TK70 <- CreateSeuratObject(counts = TK70, project = "E18",
                           min.cells = 5, min.features = 1000) 
TK70 <- subset(TK70, subset = nCount_RNA > 2000)

###merge####
TK70 <- RenameIdents(TK70, "TK70"= "E18_Rep2")
TK70 <- AddMetaData(TK70, TK70@active.ident, "Age")
head(TK70@meta.data)
Idents(TK70) <- "Age"

#check distribution of Mitochondrial genes, Ribosomal genes, and number of UMI (=nCount_RNA)
TK70[["percent.mt"]] <- PercentageFeatureSet(TK70, pattern = "^mt-")
TK70[["percent.RPS"]] <- PercentageFeatureSet(TK70, pattern = "^Rps")
TK70[["percent.RPL"]] <- PercentageFeatureSet(TK70, pattern = "^Rpl")

#filter some rough values
TK70 <- subset(TK70, subset = percent.mt <50)
TK70 <- subset(TK70, subset = percent.RPS <25)

QUALITY <- TK70@meta.data
write.csv(QUALITY, file = "CSV/E18_Rep2_QC.csv")

VlnPlot(TK70, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                           "percent.RPS","percent.RPL"), pt.size = 0)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
TK70 <- SCTransform(TK70, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
TK70 <- RunPCA(TK70, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(TK70, ndims = 50, reduction = "pca" )
#Chosen 30 to capture small changes in here

TK70 <- RunUMAP(TK70, dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(TK70, reduction = "umap", label = F,pt.size = 0.5, 
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

p <- FeaturePlot(TK70, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

TK70 <- CellCycleScoring(TK70, s.features = m.s.genes$mouseGene,
                         g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
DimPlot(TK70, group.by = "Phase")

#FeaturePlot(TK70, features = genes, cols = c("lightgrey","purple3"), pt.size = 0.5, ncol = 2, order = T )

DimPlot(TK70, reduction = "umap", label = F,pt.size = 0.5, cols =  )

####Clusters####
TK70 <- FindNeighbors(TK70, dims = 1:20)
TK70 <- FindClusters(TK70, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(TK70, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#initial markers
markers <- FindAllMarkers(TK70, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/E18_Rep2_DEG-initial.csv")

#Leave immune cells
save(TK70, file = "Robj/E18_Rep2.Robj")


#start from here again
#Merge with HyDD_v1 E18
E18 <- merge(x = TK69, y = list(TK70))
head(E18@meta.data)
E18@meta.data$percent.mt <- NULL
E18@meta.data$percent.RPS <- NULL
E18@meta.data$percent.RPL <- NULL
E18@meta.data$S.Score <- NULL
E18@meta.data$G2M.Score <- NULL
E18@meta.data$Phase <- NULL
E18@meta.data$old.ident <- NULL
E18@meta.data$SCT_snn_res.0.4 <- NULL
E18@meta.data$seurat_clusters <- NULL
save(E18, file = "Robj/E18_Merge.Robj")
load(file = "Robj/E18_Merge.Robj")
Idents(E18) <- "orig.ident"
E18[["percent.sex"]] <- PercentageFeatureSet(E18,
                                             features = c("Xist","Malat1","Tsix"))
E18[["percent.mt"]] <- PercentageFeatureSet(E18, pattern = "^mt-")
E18[["percent.RPS"]] <- PercentageFeatureSet(E18, pattern = "^Rps")
E18[["percent.RPL"]] <- PercentageFeatureSet(E18, pattern = "^Rpl")


VlnPlot(E18, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                          "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)

table(E18@meta.data$orig.ident)
Idents(E18) <- "Age"
E18<- SCTransform(E18, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA

#Run PCA
E18<- RunPCA(E18, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(E18, ndims = 50, reduction = "pca" )
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
E18 <- CellCycleScoring(E18, s.features = m.s.genes$mouseGene,
                        g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)

####harmony####
E18<- RunHarmony(E18, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
E18<- RunUMAP(E18, reduction = "harmony", dims = 1:20,
              n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(E18, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(E18, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
DimPlot(E18, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(E18, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
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
p <- FeaturePlot(E18, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = T, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

save(E18, file = "Robj/E18_Merge.Robj")

####Cluster-identity####
E18 <- FindNeighbors(E18, dims = 1:20, reduction = "harmony")
E18 <- FindClusters(E18, resolution = 2) #SCT_snn_res.0.8
DimPlot(E18, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(E18, file = "Robj/E18_Merge.Robj")

#initial markers
markers <- FindAllMarkers(E18, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/E18-Cluster-DEG-initial.csv")

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

p <- FeaturePlot(E18, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = F, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)



#Rename
#E18 what is happneing with prethal
#Where is LH
#Losing spatial pattern markers
#Huge increzse in neuropeptide/neurotransmitters
#0 VMH
#1 Tac1 SMN
#2 Junk
#3 Cck (MMN)
#4 Junk
#5 Junk? Isl1
#6 SMN
#7 Avp/Tac1
#8 Junk
#9 Cck (MMN)
#10 Junk? Isl1
#11 Junk?
#12 Glial progenitor
#13 Cck (MMN)
#14 Prethal?
#15 Tac1/Calb2
#16 Pomc/Gal
#17 Glial progenitor
#18 Glial progenitor
#19 Immune
#20 Blood
#21 Otp
#22 Prethal
#23 Glial progenitor
#24 Glial progenitor
#25 Sst/Npy/Otp
#26 Prethal?
#27 Glial progenitor
#28 Avp/Rorb/Gal (SCN)
#29 Endothelial
#30 MMN
#31 Glial progenitor
#32 Blood
#33 Npy/Agrp/Sst/Otp
#34 Junk?
#35 Glial progenitor
#36 Immune
#37 Glial progenitor
#38 Glial progenitor
#39 Immune
#40 Immune
#41 Cck (MMN)
#42 Avp/Gal/Oxt
#43 Blood
#44 Endothelial
#45 Immune


E18 <- RenameIdents(E18, "0" = "VMH", "1" = "Tac1 (SMN)", "2" = "Junk", "3" = "Cck (MMN)",
                    "4" = "Junk", "5" = "Check (ISl1/PreThal)", "6" = "SMN", "7" = "Tac1/Avp",
                    "8" = "Junk", "9" = "Cck (MMN)", "10" = "Check (ISl1/PreThal)", "11" = "Junk",
                    "12" = "Glial Pro?", "13" = "Cck (MMN)", "14" = "Check (PreThal)", 
                    "15" = "Tac1/Calb2", "16" = "Pomc/Gal", "17" = "Glial Pro?", "18" = "Glial Pro?",
                    "19" = "Immune", "20" = "Blood", "21" = "Otp", "22" = "Check (PreThal)",
                    "23" = "Glial Pro?", "24" = "Glial Pro?", "25" = "Sst/Npy/Otp", "26" = "Check (PreThal)",
                    "27" = "Glial Pro?", "28" = "Avp/Rorb/Gal (SCN?)", "29" = "Endothelial", "30" = "MMN",
                    "31" = "Glial Pro?", "32" = "Blood", "33" = "Npy/Agrp/Sst/Otp", "34" = "Check (Junk?)",
                    "35" = "Glial Pro?", "36" = "Immune", "37" = "Glial Pro?", "38" = "Glial Pro?",
                    "39" = "Immune", "40" = "Glial Pro?", "41" = "Cck (MMN)", "42" = "Avp/Gal/Oxt",
                    "43" = "Blood", "44" = "Endothelial", "45" = "Immune")
Group <- sort(levels(E18@active.ident))
E18 <- AddMetaData(E18, E18@active.ident, "Cluster_Pass1")
E18@active.ident <- factor(x = E18@active.ident, levels = Group)
E18@meta.data$Cluster_Pass1 <- factor(x = E18@meta.data$Cluster_Pass1, levels = Group)
DimPlot(E18, reduction = "umap", label = F, pt.size = 0.5)

save(E18, file = "Robj/E18_Merge.Robj")

