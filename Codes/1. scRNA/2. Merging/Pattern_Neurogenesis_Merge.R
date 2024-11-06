library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(future)
library(harmony)
set.seed(1234)
plan("multicore", workers = 10)
plan()
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram

setwd("/media/thomaskim/Data/")

####load####
load("scRNA/Robj/pattern.Robj")
load("scRNA/Robj/neurogenesis.Robj")

pattern <- UpdateSeuratObject(pattern)
neurogenesis <- UpdateSeuratObject(neurogenesis)

####Merge####
HyDD <- merge(x = pattern, y = neurogenesis)
rm(pattern)
rm(neurogenesis)
Idents(HyDD) <- "Cluster_Pass2"
#REMOVE Check (doublet #1-#4), POA_SCN?

HyDD <- subset(HyDD, idents = c("Ant_ID", "AntID_ID", "ARC (Npy_Agrp)", "ARC (Pomc)", "ARC_VMH", 
                           "ARC_VMH (Tac2_PMN)", "Astrocytes", 
                           "DMH (Npw_PMN)", "DMH-LH (Grp_Cck)", "DMH-LH (Npvf)", "Ependymal", 
                           "Isl1 cells", "LH (Hcrt_Oxt)", "LH (Pmch_Trh)", "MMN", "MMN (Cck)", 
                           "MMN (Npy_Cck)", "MMN (Nts_Tac1_Cck)", "MMN (Prog)", "Neural Pro (Ascl1)", 
                           "Neural Pro (Ascl1/Neurog2)", "Neural Pro (Neurog2)", "NPC", 
                           "NPC (Ascl1)", "NPC (Epen)", "NPC (Epen/Tan)", "NPC (Glial_Astro)", 
                           "NPC (Glial_Check)", "NPC (Glial_Oligo)", "NPC (Glial)", 
                           "NPC (gliogenic?)", "NPC (Neurog2)", "NPC (OPC)", "NSC", 
                           "NSC (Neurog2)", "NSC_NPC", "Oligodendrocytes (newly)", "OPC", 
                           "PMN", "PMN (Ghrh_Gal_Cited1)", "PMN (Ghrh_Gal_Oxt)", "PMN (Pro)", 
                           "PMN? (Gsx1)", "POA_SCN",  "PreThal", "PreThal (Pro)", 
                           "PreThal (Sp8)", "PreThal (Sp9)", "PreThal_AntID", "PreThal_ID", 
                           "PVH_SON", "PVH_SON_AntVen", "PVN_SON (Avp_Gal_Oxt)", "SCN (Rorb)", 
                           "SMN", "SMN (Calb2)", "SMN (Prog)", "SMN (Tac2)", "SST (DMH?)", 
                           "Sst_Npy", "Tanycyte", "TMN_PH (Hdc)", "Tub (ARC_VMH)", "Tub (LH-Pmch)", 
                           "Tub (LH)", "Tub (Prog)", "ZLI?"))
HyDD <- SCTransform(HyDD,
                    vars.to.regress = c("nCount_RNA","nFeature_RNA"), conserve.memory = T, ncells = 20000)
HyDD <- RunPCA(HyDD, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(HyDD, ndims = 50, reduction = "pca" )
HyDD <- RunHarmony(HyDD, group.by.vars = "orig.ident", assay.use = "SCT", plot.convergence = T)
ElbowPlot(HyDD, ndims = 50, reduction = "harmony" )

Idents(HyDD) <- "Age_Sum"
HyDD <- RunUMAP(HyDD, reduction = "harmony", dims = 1:30,
                n.neighbours = 20L, min.dist = 0.01, spread = 5)
DimPlot(HyDD, reduction = "umap", label = F, pt.size = 0.1, raster = T) + NoLegend() + NoAxes()

# Define mellow shades of blue and red
start_red = "#FF0000"   # Bright red
mid_yellow = "#FFFF00"  # Bright yellow
mid_green = "#00FF00"   # Bright green
end_purple = "#800080"  # Purple

# Create a custom gradient
colors <- colorRampPalette(c(start_red, mid_yellow, mid_green, end_purple))(9)

DimPlot(HyDD, reduction = "umap", label = T, pt.size = 0.05, raster = F, group.by = "Cluster_Pass2") + NoLegend()
DimPlot(HyDD, reduction = "umap", label = F, pt.size = 0.05, raster = F, cols = colors, group.by = "Age_Sum",) + NoLegend()
DimPlot(HyDD, reduction = "umap", label = F, pt.size = 0.05, raster = F, cols = colors, split.by = "Age_Sum",
        group.by = "Age_Sum", ncol = 3) + NoLegend()

FeaturePlot(HyDD, "Slc1a3", order = T)

save(HyDD, file = "scRNA/Robj/HyDD.Robj")


####Clean up####
HyDD <- FindNeighbors(HyDD, dims = 1:30, reduction = "harmony")
HyDD <- FindClusters(HyDD, resolution = 0.8) #SCT_snn_res.0.8
DimPlot(HyDD, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

markers <- FindAllMarkers(HyDD, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "scRNA/DEG-HyDD.csv")

#34 clusters
#31 contamnation 
#32 contamnation
FeaturePlot(HyDD, "Notch2", order = T)

HyDD <- subset(HyDD, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 
                                "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", 
                                "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "33", "34"))
HyDD <- SCTransform(HyDD,
                    vars.to.regress = c("nCount_RNA","nFeature_RNA"), conserve.memory = T, ncells = 20000)
HyDD <- RunPCA(HyDD, npcs = 50, ndims.print = NA, verbose = F)
HyDD <- RunHarmony(HyDD, group.by.vars = "Age_Sum", assay.use = "SCT", plot.convergence = T)
ElbowPlot(HyDD, ndims = 50, reduction = "harmony" )
Idents(HyDD) <- "Age_Sum"
HyDD <- RunUMAP(HyDD, reduction = "harmony", dims = 1:30,
                n.neighbours = 20L, min.dist = 0.01, spread = 5)
DimPlot(HyDD, reduction = "umap", label = F, pt.size = 0.1, raster = T) + NoLegend() + NoAxes()

# Define mellow shades of blue and red
start_red = "#FF0000"   # Bright red
mid_yellow = "#FFFF00"  # Bright yellow
mid_green = "#00FF00"   # Bright green
end_purple = "#800080"  # Purple

Mellow_Red = "#FF6666" #(Soft Coral)
Mellow_Yellow = "#FFFF66"#(Creamy Yellow)
Mellow_Green = "#66FF66" #(Mint Green)
Mellow_Purple = "#B366FF" #(Lavender)

#colors <- colorRampPalette(c(start_red, mid_yellow, mid_green, end_purple))(9)
colors <- colorRampPalette(c(Mellow_Red, Mellow_Yellow, Mellow_Green, Mellow_Purple))(9)

DimPlot(HyDD, reduction = "umap", label = T, pt.size = 0.05, raster = F, group.by = "Cluster_Pass2") + NoLegend()
DimPlot(HyDD, reduction = "umap", label = F, pt.size = 0.05, raster = F, cols = colors, group.by = "Age_Sum",) + NoLegend()
DimPlot(HyDD, reduction = "umap", label = F, pt.size = 0.05, raster = F, cols = colors, split.by = "Age_Sum",
        group.by = "Age_Sum", ncol = 3) + NoLegend()

FeaturePlot(HyDD, "Nkx2-1", order = T)

HyDD <- FindNeighbors(HyDD, dims = 1:30, reduction = "harmony")
HyDD <- FindClusters(HyDD, resolution = 0.8) #SCT_snn_res.0.8
DimPlot(HyDD, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(HyDD, file = "scRNA/Robj/HyDD.Robj")

markers <- FindAllMarkers(HyDD, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "scRNA/DEG-HyDD.csv")

Idents(HyDD) <- "Cluster_Pass2"
markers <- FindAllMarkers(HyDD, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "scRNA/DEG-HyDD_orig.cluster.csv")

####cleanup####

#17,#35 - v2 bias, #22, #26, #37 - duplicate bias
#36, #38 dissociation bias
table(HyDD@active.ident, HyDD$Age)

HyDD <- subset(HyDD, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                "11", "12", "13", "14", "15", "16", "18", "19", "20", 
                                "21", "23", "24", "25", "27", "28", "29", "30", 
                                "31", "32", "33", "34"))
HyDD <- SCTransform(HyDD,
                    vars.to.regress = c("nCount_RNA","nFeature_RNA"), conserve.memory = T, ncells = 20000)

HyDD <- RunPCA(HyDD, npcs = 50, ndims.print = NA, verbose = F)
HyDD <- RunHarmony(HyDD, group.by.vars = "orig.ident", assay.use = "SCT", plot.convergence = T)
Idents(HyDD) <- "Age_Sum"
HyDD <- RunUMAP(HyDD, reduction = "harmony", dims = 1:30,
                n.neighbours = 20L, min.dist = 0.01, spread = 5)
DimPlot(HyDD, reduction = "umap", label = F, pt.size = 0.1, raster = T) + NoLegend() + NoAxes()

Mellow_Red = "#FF6666" #(Soft Coral)
Mellow_Yellow = "#FFFF66"#(Creamy Yellow)
Mellow_Green = "#66FF66" #(Mint Green)
Mellow_Purple = "#B366FF" #(Lavender)

#colors <- colorRampPalette(c(start_red, mid_yellow, mid_green, end_purple))(9)
colors <- colorRampPalette(c(Mellow_Red, Mellow_Yellow, Mellow_Green, Mellow_Purple))(9)

DimPlot(HyDD, reduction = "umap", label = T, pt.size = 0.05, raster = F, group.by = "Cluster_Pass2") + NoLegend()
DimPlot(HyDD, reduction = "umap", label = F, pt.size = 0.05, raster = F, cols = colors, group.by = "Age_Sum",) + NoLegend()
DimPlot(HyDD, reduction = "umap", label = F, pt.size = 0.05, raster = F, cols = colors, split.by = "Age_Sum",
        group.by = "Age_Sum", ncol = 3) + NoLegend()


HyDD <- FindNeighbors(HyDD, dims = 1:30, reduction = "harmony")
HyDD <- FindClusters(HyDD, resolution = 0.8) #SCT_snn_res.0.8
DimPlot(HyDD, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(HyDD, file = "scRNA/Robj/HyDD.Robj")

markers <- FindAllMarkers(HyDD, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "scRNA/DEG-HyDD.csv")

Idents(HyDD) <- "Cluster_Pass2"
markers <- FindAllMarkers(HyDD, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "scRNA/DEG-HyDD_orig.cluster.csv")



####Identify####

#Idents(HyDD) <- "Cluster_Pass2"
Cells1 <-WhichCells(HyDD, idents = "E11")
Cells2 <-WhichCells(HyDD, idents = "E12")
Cells3 <-WhichCells(HyDD, idents = "E13")
Cells4 <-WhichCells(HyDD, idents = "E14")
Cells5 <-WhichCells(HyDD, idents = "E15")
Cells6 <-WhichCells(HyDD, idents = "E16")
Cells7 <-WhichCells(HyDD, idents = "E18")
Cells8 <-WhichCells(HyDD, idents = "P4")
Cells9 <-WhichCells(HyDD, idents = "P8")

DimPlot(HyDD, reduction = "umap", label = F, 
        pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
  ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf"))
