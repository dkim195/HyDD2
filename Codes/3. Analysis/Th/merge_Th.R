# Load required libraries
library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(future)
library(harmony)
library(stringr)
library(reshape2)
library(colorspace)
library(grid)
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram
setwd("/media/thomaskim/Data/")

####pattern####
load(file = "scRNA/Robj/pattern.Robj")

Idents(pattern) <- "Cluster_Pass2"
pattern <- RenameIdents(pattern, "Ant_ID" = "AntID_ID",
                        "NPC (gliogenic?)" = "NPC (gliogenic)",
                        "PMN? (Gsx1)" = "PMN (Gsx1)",
                        "POA/SCN?" = "POA/SCN",
                        "SST (DMH?)" = "SST (DMH)",
                        "ZLI?" = "ZLI")
identities <- Idents(pattern)
new_order <- sort(levels(identities))
Idents(pattern) <- factor(identities, levels = new_order)

VlnPlot(pattern, features = "Th")
pattern_th_cells <- subset(pattern, subset = Th > 1)
active_th_cells <- WhichCells(object = pattern, expression = Th >1, slot = 'counts')

####neurogenesis-scRNA clusters####
load(file = "scRNA/Robj/neurogenesis.Robj")
Idents(neurogenesis) <- "Cluster_Pass2"
neurogenesis <- subset(neurogenesis, idents = c("ARC (Npy_Agrp)", "ARC (Pomc)", "ARC_VMH", "ARC_VMH (Tac2_PMN)", 
                                                "Astrocytes",  "DMH-LH (Grp_Cck)", 
                                                "DMH-LH (Npvf)", "DMH (Npw_PMN)", "Ependymal", "Isl1 cells", 
                                                "LH (Hcrt_Oxt)", "LH (Pmch_Trh)", "MMN", "MMN (Cck)", "MMN (Npy_Cck)", 
                                                "MMN (Nts_Tac1_Cck)", "Neural Pro (Ascl1/Neurog2)", "Neural Pro (Neurog2)", 
                                                "NPC (Epen)", "NPC (Epen/Tan)", "NPC (Glial)", "NPC (Glial_Astro)", 
                                                "NPC (Glial_Check)", "NPC (Glial_Oligo)", "NPC (OPC)", "Oligodendrocytes (newly)", 
                                                "OPC", "PMN (Ghrh_Gal_Cited1)", "PMN (Ghrh_Gal_Oxt)", "POA_SCN", 
                                                "PreThal", "PreThal_ID", "PVN_SON (Avp_Gal_Oxt)", "SCN (Rorb)", 
                                                "SMN", "SMN (Calb2)", "SMN (Tac2)", "Sst_Npy", "Tanycyte", 
                                                "TMN_PH (Hdc)"))
neurogenesis <- RenameIdents(neurogenesis,
                             "NPC (Epen)" = "NPC (Glial)",
                             "NPC (Epen/Tan)" = "NPC (Glial)", 
                             "NPC (Glial)" = "NPC (Glial)",
                             "NPC (Glial_Astro)" = "NPC (Glial)", 
                             "NPC (Glial_Check)" = "NPC (Glial)",
                             "NPC (Glial_Oligo)" = "NPC (Glial)",
                             "NPC (OPC)" = "NPC (Glial)")
identities <- Idents(neurogenesis)
new_order <- sort(levels(identities))
Idents(neurogenesis) <- factor(identities, levels = new_order)
VlnPlot(neurogenesis, features = "Th")
neurogenesis_th_cells <- subset(neurogenesis, subset = Th > 1)
active_th_cells <- WhichCells(object = neurogenesis, expression = Th >1, slot = 'counts')

####Merge####
HyDD_TH <- merge(x = pattern_th_cells, y = neurogenesis_th_cells)
rm(pattern)
rm(neurogenesis)
rm(pattern_th_cells)
rm(neurogenesis_th_cells)

Idents(HyDD_TH) <- "Cluster_Pass2"
HyDD_TH <- SCTransform(HyDD_TH,
                    vars.to.regress = c("nCount_RNA","nFeature_RNA"))
HyDD_TH <- RunPCA(HyDD_TH, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(HyDD_TH, ndims = 50, reduction = "pca" )
HyDD_TH <- RunHarmony(HyDD_TH, group.by.vars = "orig.ident", assay.use = "SCT", plot.convergence = T)
ElbowPlot(HyDD_TH, ndims = 50, reduction = "harmony" )

Idents(HyDD_TH) <- "Age_Sum"
HyDD_TH <- RunUMAP(HyDD_TH, reduction = "harmony", dims = 1:30,
                n.neighbours = 20L, min.dist = 0.01, spread = 5)
DimPlot(HyDD_TH, reduction = "umap", label = F, pt.size = 0.1, raster = T) + NoLegend() + NoAxes()

# Define mellow shades of blue and red
start_red = "#FF0000"   # Bright red
mid_yellow = "#FFFF00"  # Bright yellow
mid_green = "#00FF00"   # Bright green
end_purple = "#800080"  # Purple

# Create a custom gradient
colors <- colorRampPalette(c(start_red, mid_yellow, mid_green, end_purple))(9)

DimPlot(HyDD_TH, reduction = "umap", label = T, pt.size = 0.05, raster = F, group.by = "Cluster_Pass2") + NoLegend()
DimPlot(HyDD_TH, reduction = "umap", label = F, pt.size = 0.05, raster = F, cols = colors, group.by = "Age_Sum",) + NoLegend()
DimPlot(HyDD_TH, reduction = "umap", label = F, pt.size = 0.05, raster = F, cols = colors, split.by = "Age_Sum",
        group.by = "Age_Sum", ncol = 3) + NoLegend()

FeaturePlot(HyDD_TH, c("Th","Hmx3","Emx2","Dlx1","Dlx2","Isl1"), order = T, pt.size = 0.1)

#clustering
HyDD_TH <- FindNeighbors(HyDD_TH, dims = 1:30, reduction = "harmony")
HyDD_TH <- FindClusters(HyDD_TH, resolution = 0.8) #SCT_snn_res.0.8
DimPlot(HyDD_TH, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
markers <- FindAllMarkers(HyDD_TH, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "DEG-HyDD_TH.csv")

table(HyDD_TH@active.ident, HyDD_TH$Cluster_Pass2) -> GROUP1
write.csv(GROUP1, file = "GROUP1.csv")

table(HyDD_TH@active.ident, HyDD_TH$Age_Sum) -> GROUP2
write.csv(GROUP2, file = "GROUP2.csv")

#clean

#0 = HMX2, HMX3, LEF1 PMN
#1 = SIX3, ESRRG, MEIS2 PRETHAL, PATTERN
#2 = GHRH, ESR1, GAL, PMN (GHRH_GAL_OXT)
#3 = ONECUT1, DRD2, PNOC, ONECUT2 PRETHAL
#4 = AVP, CALB2, CALB1
#5 = PAX6, SP9, ZIC1 PRETHAL
#6 = SEMA3C, MEIS1, GFRA2 PRETHAL
#7 = LHX6, NKX2-1, NKX2-2 PRETHAL
#8 = FOXA2, IRX1, LMX1A, SMN
#9 = JUNK
#10 = JUNK
#11 = JUNK
#12 = TRH, SIM2, FEZF1, SIM1, OTP PVH_SON
#13 = SATB2, PRLR, SIX6, ISL CELLS
#14 = NPY, SST, GPR101, SST/NPY
#15 = NHLH1, SIM1, NEUROD2 PVH_SON
#16 = NTS, TAC1, PMN
#17 = TAC1, CRHBP,  PRETHAL
#18 = NPW, DMH


HyDD_TH <- subset(HyDD_TH, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8",
                                      "12", "13", "14", "15", "16",
                                      "17", "18"))
HyDD_TH <- SCTransform(HyDD_TH,
                       vars.to.regress = c("nCount_RNA","nFeature_RNA"))
HyDD_TH <- RunPCA(HyDD_TH, npcs = 50, ndims.print = NA, verbose = F)
HyDD_TH <- RunHarmony(HyDD_TH, group.by.vars = "orig.ident", assay.use = "SCT", plot.convergence = T)
HyDD_TH <- RunUMAP(HyDD_TH, reduction = "harmony", dims = 1:30,
                   n.neighbours = 20L, min.dist = 0.01, spread = 5)
HyDD_TH <- FindNeighbors(HyDD_TH, dims = 1:30, reduction = "harmony")
HyDD_TH <- FindClusters(HyDD_TH, resolution = 0.5) #SCT_snn_res.0.8
DimPlot(HyDD_TH, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

markers <- FindAllMarkers(HyDD_TH, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "DEG-HyDD_TH.csv")

table(HyDD_TH@active.ident, HyDD_TH$Cluster_Pass2) -> GROUP1
write.csv(GROUP1, file = "GROUP1.csv")

table(HyDD_TH@active.ident, HyDD_TH$Age_Sum) -> GROUP2
write.csv(GROUP2, file = "GROUP2.csv")

#0 = ONECUT1, ONECUT2, ONECUT3, PNOC, POU6F2 - PreThal, Pattern
#1 = SIX3, ESRRG, MEIS2 - PRETHAL, Pattern
#2 = GHRH, SLC18A3, ESR1, GAL -  PMN, neurogenesis
#3 = LHX6, SP9, PAX6 - PRETHAL, Pattern
#4 = HMX2, HMX3, LEF1 - PMN, Pattern
#5 = AVP, CALB1, TAC1 - PMN, neurogenesis
#6 = FOXA2, LMX1A, N45A2 - SMN, Pattern
#7 = NTS, CITED1, TAC1, HMX2, HMX3 -  PMN, neurogenesis
#8 = SATB2, PRLR, SIX6 - ISL1, neurogenesis
#9 = SIM2, TRH, FEZF1, SIM1, OTP - PVN/SON, Pattern
#10 = GPR101, NPY, SST - CHECK, neurogenesis
#11 = NEUROD2, ROBO3, NHLH1, SIM1, OTP - PVN/SON, Pattern
#12 = NPW, NTSR2 - DMH, neurogenesis
#13 = NEUROG1, BULB1B, PRC1, EXCL

HyDD_TH <- subset(HyDD_TH, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8",
                                      "9","10","11","12"))
HyDD_TH <- SCTransform(HyDD_TH,
                       vars.to.regress = c("nCount_RNA","nFeature_RNA"))
HyDD_TH <- RunPCA(HyDD_TH, npcs = 50, ndims.print = NA, verbose = F)
HyDD_TH <- RunHarmony(HyDD_TH, group.by.vars = "orig.ident", assay.use = "SCT", plot.convergence = T)
HyDD_TH <- RunUMAP(HyDD_TH, reduction = "harmony", dims = 1:30,
                   n.neighbours = 20L, min.dist = 0.01, spread = 5)
DimPlot(HyDD_TH, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

HyDD_TH <- RenameIdents(HyDD_TH, "0" = "PreThal1", "1" = "PreThal2", "2" = "PMN_GHRH", "3" = "PreThal3",
                        "4" = "PMN", "5" = "PMN_AVP", "6" = "SMN", "7" = "PMN_NTS",
                        "8" = "ISL1+", "9" = "PVN_SON1", "10" = "GPR101_CHECK", "11" = "PVH_SON2", "12" = "DMH (NPW)")
HyDD_TH <- AddMetaData(HyDD_TH, HyDD_TH@active.ident, "TH_Type")
DimPlot(HyDD_TH, reduction = "umap", label = F, pt.size = 0.05, split.by = "TH_Type", ncol = 2)
DimPlot(HyDD_TH, reduction = "umap", label = T, pt.size = 0.05) + NoAxes() + NoLegend()



#change
#1) DMH (NPW) = NPY, AGRP - perhaps in the ARC as well
#2) GPR101_CHECK - ARC?, BARX2,Tbx3
#3) ISL1+ - Prlr, Slc6a3, Satb2
#4) PMN - HMX2, HMX3, LEF1, GSX1
#5) PMN_AVP - AVP, TAC1, GAL -> MIGHT BE PVN
#6) PMN_GHRH - GHRH, GAL, GSX1, Slc18a3 LIKELY TO BE ARC DUE TO SLC17A3
#7) PMN_NTS - Nts, Cited1, Tnr, Tac1, Hmx22, Hmx3 MORE LATERAL
#8) PreThal1 - Pou6f2, Onecut1, Meis2
#9) PreThal2 - Six3, Nr2f1, Dlx6, Reln, Sp8, Isl1
#10) PreThal3 - Lhx6, Sp9, Nkx2-2,, Arx
#11) PVH_SON2 - Nhlh1, Sim1, Otp, Lhx5, Fezf1
#12) PVN_SON1 - Trh, Fezf1, Sim1, Otp, Lhx2
#13) SMN - Foxa2, Lmx1a, Irx1, Nr4a2
HyDD_TH <- RenameIdents(HyDD_TH, "PreThal1"="TH1", "PreThal2"="TH2", "PMN_GHRH"="TH3", "PreThal3"="TH4",
                        "PMN"="TH5", "PMN_AVP"="TH6", "SMN"="TH7", "PMN_NTS"="TH8",
                        "ISL1+"="TH9", "PVN_SON1"="TH10", "GPR101_CHECK"="TH11", "PVH_SON2"="TH12", "DMH (NPW)"="TH13")
HyDD_TH <- AddMetaData(HyDD_TH, HyDD_TH@active.ident, "TH_Type2")
DimPlot(HyDD_TH, reduction = "umap", label = F, pt.size = 1) + NoAxes() + NoLegend()
DimPlot(HyDD_TH, reduction = "umap", label = T, pt.size = 0.5) + NoAxes() + NoLegend()

Idents(HyDD_TH) <- "TH_Type"
DimPlot(HyDD_TH, reduction = "umap", label = F,
        pt.size = 1) + NoAxes() + NoLegend()

DimPlot(HyDD_TH, reduction = "umap", label = T,
        pt.size = 0.5) + NoAxes() + NoLegend()


markers <- FindAllMarkers(HyDD_TH, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "Figures/Th/DEG-HyDD_TH.csv")
save(HyDD_TH, file = "scRNA/Robj/HyDD_Th.Robj")


# List of genes to plot
genes_to_plot <- c("Onecut1","Onecut2","Onecut3","Pnoc","Pou6f2","Six3","Esrrg","Meis2",
                   "Ghrh","Slc18a3","Esr1","Gal","Lhx6","Sp8","Sp8","Arx","Pax6",
                   "Hmx2","Hmx3","Lef1","Gsx1","Avp","Calb1","Tac1","Foxa2","Lmx1a",
                   "Nr4a2","Satb2","Prlr","Ar","Six6","Isl1","Npy","Gpr101","Trh","Sim1",
                   "Otp","Fezf1","Neurod2","Robo3","Nhlh1","Npw","Ntsr2",
                   "Agrp","Tbx3","Slc6a3","Slc18a3","Nts","Cited1","Tac1","Reln",
                   "Nkx2-2","Irx1","Th","Plagl1","Pbx3","Bsx","Meis1","Nkx2-4","Dlx5",
                   "Dlx6","Emx2","Mef2c", "Npy", "Agrp", "Prlr", "Pou6f2", "Meis2",
                   "Nr2f1", "Dlx6", "Reln", "Sim1", "Otp", "Lhx5","Fezf1","Trh","Sim1","Lhx2",
                   "Nkx2-2")
genes_to_plot <- unique(genes_to_plot)

# Function to create directory if it doesn't exist
create_dir_if_not_exists <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# Loop through each gene
for (gene in genes_to_plot) {
  
  # Generate the plot
  p <- FeaturePlot(HyDD_TH, features = gene, order = TRUE, pt.size = 1) + NoLegend() + NoAxes()
  
  # Construct the file name
  filename <- paste0("Figures/Th/TH_", gene, ".png")
  
  # Create directory if it doesn't exist
  create_dir_if_not_exists(dirname(filename))
  
  # Save the plot as a PNG file
  png(filename = filename, width = 900, height = 900, res = 100)
  print(p)
  dev.off()
}

