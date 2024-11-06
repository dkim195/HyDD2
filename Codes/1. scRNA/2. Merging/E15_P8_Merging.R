library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)

####loading####
#E15-E18_Merge Merge
setwd("D:/")
load(file = "Robj/E15_Merge.Robj")
load(file = "Robj/E16_Merge.Robj")
load(file = "Robj/E18_Merge.Robj")
load(file = "Robj/P4_Merge.Robj")
load(file = "Robj/P8_Merge.Robj")

neurogenesis <- merge(x = E15, y = list(E16, E18, P4, P8))

table(neurogenesis@meta.data$orig.ident)
Idents(neurogenesis) <- "Age"
table(neurogenesis@active.ident)
#order
my_levels <- c("E15_Rep1", "E15_Rep2", "E16_Rep1", "E16_Rep2",
               "E18_Rep1", "E18_Rep2",
               "P4_Rep1", "P4_Rep2", "P4_Rep3",
               "P8_Rep1", "P8_Rep2")
neurogenesis@active.ident <- factor(x = neurogenesis@active.ident, levels = my_levels)
neurogenesis@meta.data$Age <- factor(x = neurogenesis@meta.data$Age, levels = my_levels)

neurogenesis <- RenameIdents(neurogenesis, "E15_Rep1" = "E15", "E15_Rep2" = "E15",
                             "E16_Rep1" = "E16", "E16_Rep2" = "E16",
                             "E18_Rep1" = "E18", "E18_Rep2" = "E18",
                             "P4_Rep1" = "P4" , "P4_Rep2" = "P4", "P4_Rep3" = "P4",
                             "P8_Rep1" = "P8", "P8_Rep2" = "P8")
neurogenesis <- AddMetaData(neurogenesis, neurogenesis@active.ident, "Age_Sum")
my_levels <- c("E15", "E16", "E18",  "P4",  "P8")
neurogenesis@active.ident <- factor(x = neurogenesis@active.ident, levels = my_levels)
neurogenesis@meta.data$Age_Sum <- factor(x = neurogenesis@meta.data$Age_Sum, levels = my_levels)

VlnPlot(neurogenesis, features = c("nCount_RNA","nFeature_RNA",
                              "percent.mt","percent.RPS","percent.RPL"),pt.size = 0)
save(neurogenesis, file = "Robj/neurogenesis.Robj")
table(neurogenesis@active.ident)



####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
#ncells = 80000 to save some processing time
library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("D:/")
load(file = "Robj/neurogenesis.Robj")
Idents(neurogenesis) <- "Cluster_Pass1"
table(neurogenesis@active.ident)
clusters <- levels(neurogenesis@active.ident)
fix(clusters)
neurogenesis <- RenameIdents(neurogenesis, "PMM" = "PMN", 
                             "ARC (Sst/Npy)" = "ARC (Sst_Npy)", "Check (POA/SCN)" = "Check (POA_SCN)",
                              "PVH/SON?"="PVH_SON?", 
                              "PreThal/ID"="PreThal_ID", "PreThal/ID (Glia?)"= "PreThal_ID (Glia?)",
                             "Neural Pro (SMN/MMN)"="Neural Pro (SMN_MMN)", "Hcrt/Npvf"="Hcrt_Npvf",
                             "ARC/VMH"="ARC_VMH", 
                             "Pomc/Npy/Sst (Arc)"="Pomc_Npy_Sst (Arc)",  
                             "Trh/Pnoc (POA/SCN?)"="Trh_Pnoc (POA_SCN?)", "Bsx/Cck (PMN?)"="Bsx_Cck (PMN?)", 
                             "Nts/Tac1/Cck (MMN)"="Nts_Tac1_Cck (MMN)",
                             "Rora/Rorb (POA/SCN?)"="Rora_Rorb (POA_SCN?)", "Check (SMN/PVN/SON?)"="Check (SMN_PVN_SON?)", 
                             "Penk (POA/SCN?)"="Penk (POA_SCN?)",
                             "PreThal/ID (Glial?)"="PreThal_ID (Glial?)", "Avp/Gal/Oxt (PVN/SON?)"="Avp_Gal_Oxt (PVN_SON?)", 
                              "Npy/Agrp/Sst/Otp"="Npy_Agrp_Sst_Otp", "Tac1/Avp"="Tac1_Avp",
                             "Check (ISl1/PreThal)"="Check (ISl1_PreThal)", 
                              "Pomc/Gal"="Pomc_Gal", "Tac1/Calb2"="Tac1_Calb2",
                             "Avp/Rorb/Gal (SCN?)"="Avp_Rorb_Gal (SCN?)", "Sst/Npy/Otp"="Sst_Npy_Otp", 
                              "Avp/Gal/Oxt"="Avp_Gal_Oxt")
clusters <- levels(neurogenesis@active.ident)
fix(clusters)
neurogenesis <- subset(neurogenesis, idents = c("PMN", "ARC (Sst_Npy)", "Check (POA_SCN)", "PVH_SON?", "PreThal_ID", 
                                                "PreThal_ID (Glia?)", "Neural Pro (SMN_MMN)", "Hcrt_Npvf", "ARC_VMH", 
                                                "Pomc_Npy_Sst (Arc)", "Trh_Pnoc (POA_SCN?)", "Bsx_Cck (PMN?)", 
                                                "Nts_Tac1_Cck (MMN)", "Rora_Rorb (POA_SCN?)", "Check (SMN_PVN_SON?)", 
                                                "Penk (POA_SCN?)", "PreThal_ID (Glial?)", "Avp_Gal_Oxt (PVN_SON?)", 
                                                "Npy_Agrp_Sst_Otp", "Tac1_Avp", "Check (ISl1_PreThal)", "Pomc_Gal", 
                                                "Tac1_Calb2", "Avp_Rorb_Gal (SCN?)", "Sst_Npy_Otp", "Avp_Gal_Oxt", 
                                                "VMH",  "Glio Pro?", "SMN", "Check (ID?)", "MMN", "MMN (Cck)", 
                                                "Arc?", "ARC?", "Arc (Glial?)", "SMN (Glia?)", "Pomc", "Pmch", 
                                                 "Check (MMN?)", 
                                                "Glial Pro?", "Cck (MMN)", "Glial Pro? (Tanycyte?)", 
                                                "Check (Ant ID?)", "Check (SMN?)",
                                                "Tac1 (SMN)", "Gal (PMN)", 
                                                "Check (PMN?)", "Check (PreThal)", 
                                                "Otp","Glia", "Neuron" ))


neurogenesis <- SCTransform(neurogenesis,
                       vars.to.regress = c("nCount_RNA","nFeature_RNA"), ncells = 80000)
neurogenesis <- RunPCA(neurogenesis, npcs = 50, ndims.print = NA, verbose = F)
neurogenesis <- RunHarmony(neurogenesis, group.by.vars = "orig.ident", assay.use = "SCT", plot.convergence = T)
neurogenesis <- RunUMAP(neurogenesis, reduction = "harmony", dims = 1:50,
                   n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(neurogenesis, reduction = "umap", label = T,pt.size = 0.1) + NoLegend() + NoAxes()

DimPlot(neurogenesis, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age_Sum", ncol = 2) + NoLegend()

neurogenesis <- RenameIdents(neurogenesis, "PVH_SON?" = "PVH_SON_Check",  
                                             "PreThal_ID (Glia?)"="PreThal_ID (Glia_Check)", 
                                             "Trh_Pnoc (POA_SCN?)"="Trh_Pnoc (POA_SCN_Check)",
                             "Bsx_Cck (PMN?)"="Bsx_Cck (PMN_Check)", 
                                             "Rora_Rorb (POA_SCN?)"="Rora_Rorb (POA_SCN_Check)",
                             "Check (SMN_PVN_SON?)"="Check (SMN_PVN_SON_Check)", 
                                             "Penk (POA_SCN?)"="Penk (POA_SCN_Check)",
                             "PreThal_ID (Glial?)"="PreThal_ID (Glial_Check)",
                             "Avp_Gal_Oxt (PVN_SON?)"="Avp_Gal_Oxt (PVN_SON_Check)", 
                             "Avp_Rorb_Gal (SCN?)"="Avp_Rorb_Gal (SCN_Check)",
                             "Glio Pro?"="Glio Pro_Check", "Check (ID?)"="Check (ID_Check)", 
                                             "Arc?"="Arc_Check", "ARC?"="Arc_Check",
                             "Arc (Glial?)"="Arc (Glial_Check)", "SMN (Glia?)"="SMN (Glia_Check)",
                              "Check (MMN?)"="Check (MMN_Check)", "Glial Pro?"="Glial Pro_Check",
                             "Glial Pro? (Tanycyte?)"="Glial Pro (Tanycyte_Check)", 
                                             "Check (Ant ID?)"="Check (Ant ID)",
                             "Check (SMN?)"="Check (SMN)", 
                                             "Check (PMN?)"="Check (PMN)")

clusters <- levels(neurogenesis@active.ident)
auto <- function(geneid){
  if (geneid == geneid) {
    mypath <- file.path("PNG/", paste((geneid),".png",sep=""))
    png(filename =  mypath,
        width = 1200, height = 748)
    Cells1 <-WhichCells(neurogenesis, idents = geneid)
    p <-  DimPlot(neurogenesis, reduction = "umap", label = F, 
                  pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
      ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf")) 
    print(p)
    dev.off()
  }
}
lapply(clusters, auto)
dev.off()

neurogenesis <- FindNeighbors(neurogenesis, dims = 2:20, reduction = "harmony")
neurogenesis <- FindClusters(neurogenesis, resolution = 2) #SCT_snn_res.0.8
DimPlot(neurogenesis, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

save(neurogenesis, file = "Robj/neurogenesis.Robj")
markers <- FindAllMarkers(neurogenesis, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/DEG-neurogenesis.csv")

FeaturePlot(neurogenesis, "Foxg1", order = T)

DimPlot(neurogenesis, reduction = "umap", label = F,pt.size = 0.05, split.by = "Age_Sum",ncol = 2)




####Cleanup_1####
#REMOVE JUNK FIRST
#31 Junk
#35 Blood
#45 BAM
#46 Junk
neurogenesis <- subset(neurogenesis, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
                                                "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", 
                                                "23", "24", "25", "26", "27", "28", "29", "30",  "32", "33", 
                                                "34", "36", "37", "38", "39", "40", "41", "42", "43", "44", 
                                                "47", "48", "49", "50", "51", "52", "53", "54"))
options(future.globals.maxSize = 8000 * 1024^2)
neurogenesis <- SCTransform(neurogenesis,
                       vars.to.regress = c("nCount_RNA","nFeature_RNA"), verbose = T) #,  conserve.memory = T
neurogenesis <- RunPCA(neurogenesis, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(neurogenesis, ndims = 50, reduction = "pca" )

neurogenesis <- RunHarmony(neurogenesis, group.by.vars = "orig.ident", assay.use = "SCT", plot.convergence = T)
neurogenesis <- RunUMAP(neurogenesis, reduction = "harmony", dims = 1:30,
                   n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(neurogenesis, reduction = "umap", label = T,pt.size = 0.1) + NoLegend() + NoAxes()


neurogenesis <- FindNeighbors(neurogenesis, dims = 1:30, reduction = "harmony")
neurogenesis <- FindClusters(neurogenesis, resolution = 2) #SCT_snn_res.0.8
DimPlot(neurogenesis, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(neurogenesis, file = "Robj/neurogenesis.Robj")

markers <- FindAllMarkers(neurogenesis, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/DEG-neurogenesis.csv")

FeaturePlot(neurogenesis, "Foxg1", order = T)

Idents(neurogenesis) <- "Cluster_Pass1a"
neurogenesis <- AddMetaData(neurogenesis, neurogenesis@active.ident, "Cluster_Pass1a")
clusters <- levels(neurogenesis$Cluster_Pass1a)
auto <- function(geneid){
  if (geneid == geneid) {
    mypath <- file.path("PNG/", paste((geneid),".png",sep=""))
    png(filename =  mypath,
        width = 1200, height = 748)
    Cells1 <-WhichCells(neurogenesis, idents = geneid)
    p <-  DimPlot(neurogenesis, reduction = "umap", label = F, 
                  pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
      ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf")) 
    print(p)
    dev.off()
  }
}
lapply(clusters, auto)
dev.off()


####Cleanup 2####
#Remove Junk
#0 Junk
#17 Ribosomal
#30 Junk
#48 Junk
#52 Doublet
#53 Doublet
#56 Doublet

cluster <- levels(neurogenesis@active.ident)
fix(cluster)
neurogenesis <- AddMetaData(neurogenesis, neurogenesis@active.ident, "Cluster_Pass1a")
options(future.globals.maxSize = 8000 * 1024^2)
neurogenesis <- subset(neurogenesis, idents = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
                                                  "12", "13", "14", "15", "16",  "18", "19", "20", "21", "22", 
                                                  "23", "24", "25", "26", "27", "28", "29", "31", "32", "33", 
                                                  "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", 
                                                  "45", "46", "47",  "49", "50", "51",  "54", "55","57"))
neurogenesis <- SCTransform(neurogenesis,
                            vars.to.regress = c("nCount_RNA","nFeature_RNA"), verbose = T) #,  conserve.memory = T
neurogenesis <- RunPCA(neurogenesis, npcs = 50, ndims.print = NA, verbose = F)
neurogenesis <- RunHarmony(neurogenesis, group.by.vars = "orig.ident", assay.use = "SCT", plot.convergence = T)
neurogenesis <- RunUMAP(neurogenesis, reduction = "harmony", dims = 1:30,
                        n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(neurogenesis, reduction = "umap", label = T,pt.size = 0.1) + NoLegend() + NoAxes()


neurogenesis <- FindNeighbors(neurogenesis, dims = 1:30, reduction = "harmony")
neurogenesis <- FindClusters(neurogenesis, resolution = 2) #SCT_snn_res.0.8
DimPlot(neurogenesis, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(neurogenesis, file = "Robj/neurogenesis.Robj")

markers <- FindAllMarkers(neurogenesis, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.05, verbose = T)
write.csv(markers, file = "CSV/DEG-neurogenesis.csv")

Idents(neurogenesis) <- "Cluster_Pass1a"
neurogenesis <- AddMetaData(neurogenesis, neurogenesis@active.ident, "Cluster_Pass1a")
clusters <- levels(neurogenesis$Cluster_Pass1a)
auto <- function(geneid){
  if (geneid == geneid) {
    mypath <- file.path("PNG/", paste((geneid),".png",sep=""))
    png(filename =  mypath,
        width = 1200, height = 748)
    Cells1 <-WhichCells(neurogenesis, idents = geneid)
    p <-  DimPlot(neurogenesis, reduction = "umap", label = F, 
                  pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
      ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf")) 
    print(p)
    dev.off()
  }
}
lapply(clusters, auto)
dev.off()


####Continue###

#0 OPC
#1 MMN (CCK)
#2 NPC (Epen/Tan)
#3 ARC_VMH
#4 NPC (Astro/glial)
#5 NPC (Astro/glial)
#6 POA_SCN
#7 NPC (Tanycyte)
#8 NPC (Astro/glial) more oligo
#9 SMN (Calb2)
#10 Isl1, Zdbf2
#11 PreThal_ID
#12 NPC (Astro/glial) more astro
#13 MMN
#14 SMN
#15 Pmch_Trh (LH)
#16 Neural Pro (Neurog2)
#17 MMN (Cck)
#18 Ependymal
#19 Nts_Tac1_Cck (MMN)
#20 PreThal
#21 Astro
#22 Ependymal
#23 Ghrh_Gal (PMN - express Oxt)
#24 Neural Pro (Ascl1/Neurog2)
#25 NPC (Glial-check)
#26 Oligo
#27 Npvf (DMH-LH)
#28 NPC (Glial-check)
#29 NPC (Glial-check)
#30 SMN (Tac2)
#31 MMN (CCK)
#32 Ghrh_Gal (Pmn, Cited1)
#33 Sst_Npy
#34 Avp_Gal_Oxt (PVN_SON)
#35 Sst_Npy
#36 Grp_Cck (LH)
#37 MMN (CCK)
#38 NPC (Epen)
#39 Hcrt_Oxt (LH)
#40 Npy_Agrp (ARC)
#41 NPC(Oligo-pre mylin)
#42 check (doublet?)
#43 Pomc (ARC)
#44 Npw (DMH, pmh derived))
#45 Oligo
#46 NPC (Astro/glial)
#47 Ghrh_Gal (Pmn, why Oxtocin)
#48 SCN (rORB)
#49 Pmch_Trh (LH)
#50 OPC
#51 check (doublet?) 
#52 OPC
#53 Tac2 (arc_vmh, orignated from PMN?)
#54 Npy_Cck (MMN)
#55 Hdc (TMN_PH)

clusters <- levels(neurogenesis@active.ident)
fix(clusters)
Idents(neurogenesis) <- "SCT_snn_res.2"
neurogenesis <- RenameIdents(neurogenesis, "0" = "OPC", "1" = "MMN (Cck)", "2" = "NPC (Epen/Tan)",
                             "3" = "ARC_VMH", "4" = "NPC (Glial)", "5" = "NPC (Glial)", "6" = "POA_SCN",
                             "7" = "Tanycyte", "8" = "NPC (Glial_Oligo)",
                             "9" = "SMN (Calb2)", "10" = "Isl1 cells",
                             "11" = "PreThal_ID", 
                             "12" = "NPC (Glial_Astro)", "13" = "MMN",
                             "14" = "SMN", "15" = "LH (Pmch_Trh)",
                             "16" = "Neural Pro (Neurog2)", "17" = "MMN (Cck)",
                             "18" = "Ependymal", "19" = "MMN (Nts_Tac1_Cck)",
                             "20" = "PreThal", "21" = "Astrocytes", "22" = "Ependymal", 
                             "23" = "PMN (Ghrh_Gal_Oxt)", "24" = "Neural Pro (Ascl1/Neurog2)",
                             "25" = "NPC (Glial_Check)",
                             "26" = "Oligodendrocytes (newly)", "27" = "DMH-LH (Npvf)", "28" = "NPC (Glial_Check)",
                             "29" = "NPC (Glial_Check)", "30" = "SMN (Tac2)",
                             "31" = "MMN (Cck)", "32" = "PMN (Ghrh_Gal_Cited1)", "33" = "Sst_Npy", 
                             "34" = "PVN_SON (Avp_Gal_Oxt)", "35" = "Sst_Npy",
                             "36" = "DMH-LH (Grp_Cck)", "37" = "MMN (Cck)",
                             "38" = "NPC (Epen)", "39" = "LH (Hcrt_Oxt)", "40" = "ARC (Npy_Agrp)",
                             "41" = "NPC (OPC)", "42" = "check (doublet #1)", "43" = "ARC (Pomc)", "44" = "DMH (Npw_PMN)", 
                             "45" = "Oligodendrocytes (newly)", "46" = "NPC (Glial)",
                             "47" = "PMN (Ghrh_Gal_Oxt)", "48" = "SCN (Rorb)",
                             "49" = "LH (Pmch_Trh)", "50" = "check (doublet #2)", "51" = "check (doublet #3)",
                             "52" = "check (doublet #4)", "53" = "ARC_VMH (Tac2_PMN)",
                             "54" = "MMN (Npy_Cck)", "55" = "TMN_PH (Hdc)" )
Group <- sort(levels(neurogenesis@active.ident))
neurogenesis <- AddMetaData(neurogenesis, neurogenesis@active.ident, "Cluster_Pass2")
neurogenesis@active.ident <- factor(x = neurogenesis@active.ident, levels = Group)
neurogenesis@meta.data$Cluster_Pass2 <- factor(x = neurogenesis@meta.data$Cluster_Pass2, levels = Group)
DimPlot(neurogenesis, reduction = "umap", label = T,pt.size = 0.05)

save(neurogenesis, file = "Robj/neurogenesis.Robj")

markers <- FindAllMarkers(neurogenesis, test.use = "wilcox",
                          logfc.threshold = 0.3,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/DEG-neurogenesis.csv")

#####Check!!#####
#cluster3 = ARC_VMH might be LH #39 merge with #3 vmh and lh 
#gaba and glut origin
#check lh origin
#where is LH GABA?
#read alexe jackson paper - otp ddc th cbln1, slc18a2

neurogenesis <- subset(neurogenesis, idents = c())
neurogenesis <- SCTransform(neurogenesis,
                       vars.to.regress = c("nCount_RNA","nFeature_RNA"))
neurogenesis <- RunPCA(neurogenesis, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(neurogenesis, ndims = 50, reduction = "pca" )

neurogenesis <- RunHarmony(neurogenesis, group.by.vars = "orig.ident",
                      assay.use = "SCT", plot.convergence = T)

neurogenesis <- RunUMAP(neurogenesis, reduction = "harmony", dims = 1:18,
                   n.neighbours = 10L, min.dist = 0.01, spread = 3)
DimPlot(neurogenesis, reduction = "umap",
        label = T, pt.size = 0.1) + NoLegend() + NoAxes()

save(neurogenesis, file = "Robj/neurogenesis.Robj")


markers <- FindAllMarkers(neurogenesis, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/DEG-neurogenesis.csv")








####CHECK####

library(viridis)
library(viridisLite)
library(ggplot2)
Grid <- neurogenesis@reductions$umap@cell.embeddings
Grid <- data.frame(Grid)
Grid$orig.ident <- neurogenesis$Age_Sum
Grid_ds <- Grid[,1:2]
ggplot(Grid, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(data=Grid_ds, size=0.1, alpha=0.1, color="white") +
  scale_fill_viridis(option="A", name = "Density")+
  facet_grid(~orig.ident) +
  stat_density_2d(geom="raster",aes(fill=stat(ndensity)), contour=F) +
  theme_classic() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 16),
        strip.text.x = element_text(size = 24)) +
  coord_fixed()



library(dplyr)
library(reshape2)
library(Rmisc)
Idents(neurogenesis) <- "SCT_snn_res.2"
meta.data <- neurogenesis@meta.data
graph <- 100*prop.table(table(Idents(neurogenesis), neurogenesis$Age), margin = 1) #margin 1= row, 2 = column
dat <- melt(graph)
ggplot(data = dat, aes(x=factor(Var1), y=value, fill=Var2)) + 
  geom_bar(fun.y = "mean", colour = "black", stat = "summary", width = 1) +
  ylab("Cluster proportion") + xlab("Cluster identity") +
  ggtitle("Proportion of cell types") +
  guides(fill=guide_legend(title="Group"))  +
  theme(axis.text.x =  element_text(size=10, angle = 75, hjust =1),
        axis.text.y = element_text(size=10), axis.title.y = element_text(size=15),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) 




Cells1 <-WhichCells(neurogenesis, idents = "24")
DimPlot(neurogenesis, reduction = "umap", label = F, 
        pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
  ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf")) 


genes <- c("Pmch","Lhx9","Hcrt","Igf1","Cartpt","Npvf","C1ql3","Rorb","Nkx2-2")
genes <- c("Lhx1","Lhx8","Onecut2","Onecut3","Ddc","Th","Pnoc")
genes <- c("Six3","Otp","Isl1","Islr2","Onecut2","Meis2","Dlx1","Dlx2","Arx","Npy","Pax6","Pou6f2")
p <- FeaturePlot(neurogenesis, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = F, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)
