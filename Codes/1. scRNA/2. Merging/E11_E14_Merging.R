library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)

####loading####
#E11-E13_Merge Merge
setwd("D:/")
load(file = "Robj/E11_Merge.Robj")
load(file = "Robj/E12_Merge.Robj")
load(file = "Robj/E13_Merge.Robj")
load(file = "Robj/E14_Merge.Robj")

pattern <- merge(x = E11, y = list(E12, E13, E14))

table(pattern@meta.data$orig.ident)
Idents(pattern) <- "Age"
table(pattern@active.ident)
#order
my_levels <- c("E11", "E11_Rep", "E12", "E12_Rep1", "E12_Rep2",
               "E13", "E13_Rep1", "E13_Rep2", "E14_Rep1", "E14_Rep2")
pattern@active.ident <- factor(x = pattern@active.ident, levels = my_levels)
pattern@meta.data$Age <- factor(x = pattern@meta.data$Age, levels = my_levels)

VlnPlot(pattern, features = c("nCount_RNA","nFeature_RNA",
                            "percent.mt","percent.RPS","percent.RPL"),
        pt.size = 0)
table(pattern@active.ident)

####processing####
#New normalization and scale data function
#assay SCT scale.data (pearson residual for PCA), count (corrected UMI for visualisation from scale.data),
#log1p DATA for differentiation
#filtered data in assay RNA cont
#ncells = 80000 to save some processing time
Idents(pattern) <- "Cluster_Pass1"
pattern <- RenameIdents(pattern, "Check (Telen RP/ZLI)" = "Check (Telen RP_ZLI)",
                        "Foxg1 (NPC/Thalamic eminence?)" = "Foxg1 (NPC_Thalamic eminence)",
                        "Foxg1/PVH/SON" = "Foxg1_PVH_SON", "Check (Meis2/Sst)" = "Check (Meis2_Sst)",
                        "Check (MMN/TT/ID)" = "Check (MMN_TT_ID)",
                        "Neural Pro (Neurog2/SMN/MMN)" = "Neural Pro (Neurog2_SMN_MMN)",
                        "Neural Pro (SMN/MMN)" = "Neural Pro (SMN_MMN)", "G2M NPC/NPC" = "G2M NPC_NPC",
                        "PVH/SON" = "PVH_SON", "Prethaalmus" = "Prethalamus",
                        "Tuberal (ARC/VMH)" = "Tuberal (ARC_VMH)", "Tuberal (PMN-LH)" = "Tuberal (PMN_LH)",
                        "SMN/MMN (TT?)" = "SMN_MMN (TT)", "SMN/MMN" = "SMN_MMN",
                        "Check (Junk?)" = "Check (Junk)", "ID?" = "ID",
                        "Prethalamus (Ant ID?)" = "Prethalamus (Ant ID)", "Ant ID?" = "Ant ID",
                        "Check (MGE?)" = "Check (MGE)", "Check (DMH?)" = "Check (DMH)",
                        "Check (Pax7 TT?)" = "Check (Pax7 TT)", "Check (PMN?)" = "Check (PMN)")
clusters <- levels(pattern@active.ident)
fix(clusters)
#remove junk thalamic eminence
pattern <- subset(pattern, idents = c("Foxg1_PVH_SON", 
                                      "Check (Meis2_Sst)", "Check (MMN_TT_ID)", "Neural Pro (Neurog2_SMN_MMN)", 
                                      "Neural Pro (SMN_MMN)", "G2M NPC_NPC", "PVH_SON", "Prethalamus", 
                                      "Tuberal (ARC_VMH)", "Tuberal (PMN_LH)", "SMN_MMN (TT)", "SMN_MMN", 
                                      "ID", "Prethalamus (Ant ID)", "Ant ID",  
                                      "Check (DMH)", "Check (Pax7 TT)", "Check (PMN)",
                                        "Neural Pro (Neurog1)", 
                                      "G2M NPC", "Neural Pro (Ascl1)", "NPC", "NPC (Tuberal)", 
                                      "Neural Pro", "NPC (Foxg1)", "MMN", "SMN", "Tuberal (PMN)", 
                                      "Neural Pro (Neurog2)", "Prethalamus (ID)",
                                      "PMN", "Check (MMN)", "Tuberal (LH)"))
save(pattern, file = "Robj/pattern.Robj")

#
library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("D:/")
load(file = "Robj/pattern.Robj")
pattern <- SCTransform(pattern,
                       vars.to.regress = c("nCount_RNA","nFeature_RNA"))
#Run PCA
pattern <- RunPCA(pattern, npcs = 50, ndims.print = NA, verbose = F)
pattern <- RunHarmony(pattern, group.by.vars = "orig.ident", assay.use = "SCT", plot.convergence = T)
ElbowPlot(pattern, ndims = 50, reduction = "pca" )
pattern <- RunUMAP(pattern, reduction = "harmony", dims = 1:50,
                 n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(pattern, reduction = "umap", label = T,pt.size = 0.1) + NoLegend() + NoAxes()
save(pattern, file = "Robj/pattern.Robj")

DimPlot(pattern, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(pattern, reduction = "umap", label = F,pt.size = 0.1, group.by = "Cluster_Pass1")
DimPlot(pattern, reduction = "umap", label = F,
        pt.size = 0.05, group.by = "Cluster_Pass1") + NoLegend()

####CLEANUP 0####
#remove Foxg1 (including Foxg1_PVH_SON = BNST)
library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("D:/")
load(file = "Robj/pattern.Robj")
pattern <- AddMetaData(pattern, pattern@active.ident, "Cluster_Pass1a")
pattern <- subset(pattern, idents = c("Check (Meis2_Sst)", "Check (MMN_TT_ID)",
                                      "Neural Pro (Neurog2_SMN_MMN)", 
                                      "Neural Pro (SMN_MMN)", "G2M NPC_NPC",
                                      "PVH_SON", "Prethalamus", 
                                      "Tuberal (ARC_VMH)", "Tuberal (PMN_LH)", "SMN_MMN (TT)", "SMN_MMN", 
                                      "ID", "Prethalamus (Ant ID)", "Ant ID",  
                                      "Check (DMH)", "Check (Pax7 TT)", "Check (PMN)",
                                      "Neural Pro (Neurog1)", 
                                      "G2M NPC", "Neural Pro (Ascl1)", "NPC", "NPC (Tuberal)", 
                                      "Neural Pro", "MMN", "SMN", "Tuberal (PMN)", 
                                      "Neural Pro (Neurog2)", "Prethalamus (ID)",
                                      "PMN", "Check (MMN)", "Tuberal (LH)"))
pattern <- SCTransform(pattern,
                       vars.to.regress = c("nCount_RNA","nFeature_RNA"))
pattern <- RunPCA(pattern, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(pattern, ndims = 50, reduction = "pca" )
pattern <- RunHarmony(pattern, group.by.vars = "orig.ident", assay.use = "SCT", plot.convergence = T)

pattern <- RunUMAP(pattern, reduction = "harmony", dims = 1:20,
                   n.neighbours = 10L, min.dist = 0.001, spread = 3)
DimPlot(pattern, reduction = "umap", label = T,pt.size = 0.1) + NoLegend() + NoAxes()


#20- spread 3, n 5 > 12 (too clump)

clusters <- levels(pattern@active.ident)
auto <- function(geneid){
  if (geneid == geneid) {
    mypath <- file.path("PNG/", paste((geneid),".png",sep=""))
    png(filename =  mypath,
        width = 1200, height = 748)
    Cells1 <-WhichCells(pattern, idents = geneid)
    p <-  DimPlot(pattern, reduction = "umap", label = F, 
                  pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
      ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf")) 
    print(p)
    dev.off()
  }
}
lapply(clusters, auto)
dev.off()

pattern <- FindNeighbors(pattern, dims = 1:20, reduction = "harmony")
pattern <- FindClusters(pattern, resolution = 0.8) #SCT_snn_res.0.8
DimPlot(pattern, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(pattern, file = "Robj/pattern.Robj")
markers <- FindAllMarkers(pattern, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "CSV/DEG-pattern.csv")

FeaturePlot(pattern, "Foxg1", order = T)

Idents(pattern) <- "Age"
pattern <- RenameIdents(pattern, "E11_Rep" = "E11", "E12_Rep1" = "E12",
                        "E12_Rep2" = "E12", "E13_Rep1" = "E13", "E13_Rep2" = "E13",
                        "E14_Rep1" = "E14", "E14_Rep2" = "E14")
pattern <- AddMetaData(pattern, pattern@active.ident, "Age_Sum")
DimPlot(pattern, reduction = "umap", label = F,pt.size = 0.05, split.by = "Age_Sum",ncol = 2)

####Cleanup_1####
#0 SMN
#1 mmn
#2 NPC
#3 Prethal sp8
#4 PVH/SON
#5 neural pro neurog2
#6 NPC Foxg1
#7 PMN
#8 NPC -glial OR ZLI
#9 NPC Foxg1
#10 Ant ID/ID
#11 Tuberal (LH) + TUBRAL (ARC/VMH)
#12 neurla pro ascl1
#13 Prethal sp9
#14 NPC
#15 NPC - junk?
#16 NPC
#17 PVH/SON  cartpt and POA
#18 MMN Pro
#19 NPC  Foxg1
#20 PMN + TUBRAL (ARC/VMH)
#21 DMH? SST?
#22 SMN midbrain
#23 junk
#24 PMN Gsx1
#25 LH?
#26 junk
#EXCL 15, 23, 26, 22
Idents(pattern) <- "SCT_snn_res.0.8"
pattern <- subset(pattern, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
                                        "12", "13", "14","16", "17", "18", "19", "20", "21",
                                        "24", "25"))
pattern <- SCTransform(pattern,
                       vars.to.regress = c("nCount_RNA","nFeature_RNA"))
pattern <- RunPCA(pattern, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(pattern, ndims = 50, reduction = "pca" )
pattern <- RunHarmony(pattern, group.by.vars = "orig.ident", assay.use = "SCT", plot.convergence = T)

pattern <- RunUMAP(pattern, reduction = "harmony", dims = 1:20,
                   n.neighbours = 10L, min.dist = 0.01, spread = 3)
DimPlot(pattern, reduction = "umap", label = T,pt.size = 0.1) + NoLegend() + NoAxes()


Idents(pattern) <- "Cluster_Pass1a"
pattern <- AddMetaData(pattern, pattern@active.ident, "Cluster_Pass1a")
clusters <- levels(pattern$Cluster_Pass1a)
auto <- function(geneid){
  if (geneid == geneid) {
    mypath <- file.path("PNG/", paste((geneid),".png",sep=""))
    png(filename =  mypath,
        width = 1200, height = 748)
    Cells1 <-WhichCells(pattern, idents = geneid)
    p <-  DimPlot(pattern, reduction = "umap", label = F, 
                  pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
      ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf")) 
    print(p)
    dev.off()
  }
}
lapply(clusters, auto)
dev.off()

pattern <- FindNeighbors(pattern, dims = 1:20, reduction = "harmony")
pattern <- FindClusters(pattern, resolution = 2) #SCT_snn_res.0.8
DimPlot(pattern, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
save(pattern, file = "Robj/pattern.Robj")

markers <- FindAllMarkers(pattern, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/DEG-pattern.csv")

FeaturePlot(pattern, "Foxg1", order = T)


Idents(pattern) <- "Age"
pattern <- RenameIdents(pattern, "E11_Rep" = "E11", "E12_Rep1" = "E12",
                        "E12_Rep2" = "E12", "E13_Rep1" = "E13", "E13_Rep2" = "E13",
                        "E14_Rep1" = "E14", "E14_Rep2" = "E14")
pattern <- AddMetaData(pattern, pattern@active.ident, "Age_Sum")
DimPlot(pattern, reduction = "umap", label = F,pt.size = 0.05, split.by = "Age_Sum",ncol = 2)


####Cleanup 2####
#Remove Junk
#0 Prethal - sp8
#1 PVH_SON
#2 SMN
#3 NPC Type 1
#4 MMN
#5 SMN - pROG
#6 Prethal Sp9
#7 Tuberal (ARC/VMH)
#8 AntID/ID
#9 NPC (gliogenic)?
#10 NPC
#11 NPC
#12 PMN
#13 G2M NPC  
#14 Neurog2
#15 G2M NPC  
#16 MMN pro
#17 Ascl1
#18 PVH_SON Cartpt + POA (anteroventral: "Mab21l2","Six3","Six6","Zic2","Vax1","Arx","Lhx1")
#19 ZLI?
#20 MMN
#21 SMN
#22 Junk
#23 Tub Pro
#24 G2M NPC/NPC
#25 PMN Pro
#26  Foxg1/NPC Neurog2
#27  AntID
#28 Prethal pro
#29 pvh/son
#30 G2M NPC
#31 SST DMH
#32 Prethal/antid
#33 NPC/Ascl
#34 G2M NPC
#35 SST DMH
#36 NPC (gliogenic)?
#37 PMN (Gsx1)
#38 Tuberal (LH)
#39 Tuberal (LH - PMCH)
#40 G2M NPC 
#41 SCN?
#42 Junk

Idents(pattern) <- "SCT_snn_res.2"
pattern <- RenameIdents(pattern, "0" = "PreThal (Sp8)", "1" = "PVH_SON",
                        "2" = "SMN", "3" = "NPC (Neurog2)", "4" = "MMN",
                        "5" = "SMN (Prog)", "6" = "PreThal (Sp9)",
                        "7" = "Tub (ARC_VMH)", "8" = "AntID_ID",
                        "9" = "NPC (gliogenic?)", "10" = "NPC", "11" = "NPC", 
                        "12" = "PMN", "13" = "G2M NPC", "14" = "Neural Pro (Neurog2)",
                        "15" = "G2M NPC", "16" = "MMN (Prog)", "17" = "Neural Pro (Ascl1)",
                        "18" = "PVH_SON_AntVen", "19" = "ZLI?",
                        "20" = "MMN", "21" = "SMN", "22" = "Junk", 
                        "23" = "Tub (Prog)", "24" = "G2M NPC_NPC", "25" = "PMN (Pro)",
                        "26" = "G2M NPC (Neurog2)", "27" = "Ant_ID", "28" = "PreThal (Pro)",
                        "29" = "PVH_SON", "30" = "G2M NPC", "31" = "SST (DMH?)",
                        "32" = "PreThal_AntID", "33" = "NPC (Ascl1)", 
                        "34" = "G2M NPC", "35" = "SST (DMH?)",
                        "36" = "NPC (gliogenic?)", "37" = "PMN? (Gsx1)",
                        "38" = "Tub (LH)", "39" = "Tub (LH-Pmch)",
                        "40" = "G2M NPC", "41" = "POA/SCN?", "42" = "Junk")
Group <- sort(levels(pattern@active.ident))
pattern <- AddMetaData(pattern, pattern@active.ident, "Cluster_Pass2")
pattern@active.ident <- factor(x = pattern@active.ident, levels = Group)
pattern@meta.data$Cluster_Pass2 <- factor(x = pattern@meta.data$Cluster_Pass2, levels = Group)
DimPlot(pattern, reduction = "umap", label = T,pt.size = 0.05)

pattern <- subset(pattern, idents = c("Ant_ID", "AntID_ID", "MMN", "MMN (Prog)", "Neural Pro (Ascl1)", 
                                        "Neural Pro (Neurog2)", "NPC", "NPC (Ascl1)", "NPC (gliogenic?)", 
                                        "NPC (Neurog2)", "G2M NPC", "G2M NPC (Neurog2)", "G2M NPC_NPC", "PMN", "PMN (Pro)", 
                                        "PMN? (Gsx1)", "POA/SCN?", "PreThal (Pro)", "PreThal (Sp8)", 
                                        "PreThal (Sp9)", "PreThal_AntID", "PVH_SON", "PVH_SON_AntVen", 
                                        "SMN", "SMN (Prog)", "SST (DMH?)", "Tub (ARC_VMH)", "Tub (LH-Pmch)", 
                                        "Tub (LH)", "Tub (Prog)", "ZLI?"))
pattern <- SCTransform(pattern,
                       vars.to.regress = c("nCount_RNA","nFeature_RNA"))
pattern <- RunPCA(pattern, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(pattern, ndims = 50, reduction = "pca" )

pattern <- RunHarmony(pattern, group.by.vars = "orig.ident",
                      assay.use = "SCT", plot.convergence = T)

pattern <- RunUMAP(pattern, reduction = "harmony", dims = 1:18,
                   n.neighbours = 10L, min.dist = 0.01, spread = 3)
DimPlot(pattern, reduction = "umap",
        label = T, pt.size = 0.1) + NoLegend() + NoAxes()

save(pattern, file = "Robj/pattern.Robj")


markers <- FindAllMarkers(pattern, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "CSV/DEG-pattern.csv")

#####START FROM HERE####
#Leave ZLI for now
#PreThal_AntID = PreThal_DorsalID, is ZI hybrid or reticular thalamus A13 and LhX6 ID
#AntID = SCN
#POA/SCN? = POA
#NPC Gliogenic = Gliogenic + Radial glia
#Early vs late NPC?
#

####CHECK####

library(viridis)
library(viridisLite)
library(ggplot2)
Grid <- pattern@reductions$umap@cell.embeddings
Grid <- data.frame(Grid)
Grid$orig.ident <- pattern$Age_Sum
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
Idents(pattern) <- "SCT_snn_res.2"
meta.data <- pattern@meta.data
graph <- 100*prop.table(table(Idents(pattern), pattern$Age), margin = 1) #margin 1= row, 2 = column
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




Cells1 <-WhichCells(pattern, idents = "24")
DimPlot(pattern, reduction = "umap", label = F, 
        pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
  ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf")) 


genes <- c("Pmch","Lhx9","Hcrt","Igf1","Cartpt","Npvf","C1ql3","Rorb","Nkx2-2")
genes <- c("Lhx1","Lhx8","Onecut2","Onecut3","Ddc","Th","Pnoc")
genes <- c("Six3","Otp","Isl1","Islr2","Onecut2","Meis2","Dlx1","Dlx2","Arx","Npy","Pax6","Pou6f2")
p <- FeaturePlot(pattern, features =  genes, 
                 cols = c("lightgrey","purple3"), pt.size = 0.5, order = F, combine = F) 
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)



