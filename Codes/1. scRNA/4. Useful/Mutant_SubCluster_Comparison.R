library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("D:/")
#load(file = "Robj/Isl1_Final.Robj")
Idents(Mutant) <- "Cluster_Pass2"
levels(Mutant@active.ident)
Mutant_Sub <- subset(Mutant, idents = c("Avp_Gal_Oxt (PVN_SON)"))
Mutant_Sub <- subset(Mutant, idents = c("Check (POA_SCN)"))
Mutant_Sub <- subset(Mutant, idents = c("Rora_Rorb (POA_SCN)"))
Mutant_Sub <- subset(Mutant, idents = c("NPC (Tanycyte)"))
Mutant_Sub <- subset(Mutant, idents = c("PreThal"))

Mutant_Sub <- SCTransform(Mutant_Sub, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
Mutant_Sub<- RunPCA(Mutant_Sub, npcs = 50, ndims.print = NA, verbose = F)
Mutant_Sub<- RunHarmony(Mutant_Sub, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
Mutant_Sub<- RunUMAP(Mutant_Sub, reduction = "harmony", dims = 1:22,
                 n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(Mutant_Sub, split.by = "Genotype")
Idents(Mutant_Sub) <- "Genotype"
markers <- FindAllMarkers(Mutant_Sub,
                       test.use = "LR", 
                       latent.vars = c("nCount_RNA","nFeature_RNA"),
                       logfc.threshold = 0.1,
                       min.pct = 0.05, verbose = T)
write.csv(markers, "Pre_Markers.csv")
View(markers)


FeaturePlot(E12_Atlas, "Nhp2l1")
FeaturePlot(Mutant_Sub, c("Pnoc","Sp8","Sp9","Dlx1","Gad2","Gad1"), split.by = "Genotype")
VlnPlot(Mutant_Sub, c("Pnoc","Sp8","Sp9","Dlx1","Gad2","Gad1"), group.by = "Genotype")

Cells1 <-WhichCells(Mutant_Sub, idents = "Prethalamus (Sp8)")
Cells1 <-WhichCells(Mutant_Sub, idents = "Prethalamus (Sp9)")
DimPlot(Mutant_Sub, reduction = "umap", label = F, 
        pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
  ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf"))
