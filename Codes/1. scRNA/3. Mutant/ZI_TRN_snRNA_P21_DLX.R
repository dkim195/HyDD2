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
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram
setwd("/media/thomaskim/")

#setwd("/media/thomaskim/Data/")
#snRNA P21, ZI/TRN

####Loading####
#Foxd1;DLX Het
TKD7<- Read10X_h5("HyDD2-HD2/scRNA/TKD7/raw_feature_bc_matrix.h5")
colnames(TKD7) = paste0("TKD7_", colnames(TKD7))
TKD7<- CreateSeuratObject(counts = TKD7, project = "DLX",
                          min.cells = 5, min.features = 1000) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TKD7<- subset(TKD7, subset = nCount_RNA > 2000)  
TKD7<- RenameIdents(TKD7, "TKD7"= "Ctrl1")
TKD7<- AddMetaData(TKD7, TKD7@active.ident, "Genotype")
Idents(TKD7) <- "Genotype"
TKD7[["percent.mt"]] <- PercentageFeatureSet(TKD7, pattern = "^mt-")
TKD7[["percent.RPS"]] <- PercentageFeatureSet(TKD7, pattern = "^Rps")
TKD7[["percent.RPL"]] <- PercentageFeatureSet(TKD7, pattern = "^Rpl")
TKD7<- subset(TKD7, subset = percent.mt <50)
TKD7<- subset(TKD7, subset = percent.RPS <25)

#Foxd1;DLX Het
TKF2 <- Read10X_h5("HyDD2-HD2/scRNA/TKF2/raw_feature_bc_matrix.h5")
colnames(TKF2) = paste0("TKF2_", colnames(TKF2))
TKF2 <- CreateSeuratObject(counts = TKF2, project = "DLX",
                           min.cells = 5, min.features = 1000) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TKF2 <- subset(TKF2, subset = nCount_RNA > 2000)  
TKF2 <- RenameIdents(TKF2, "TKF2"= "Ctrl2")
TKF2 <- AddMetaData(TKF2, TKF2@active.ident, "Genotype")
Idents(TKF2) <- "Genotype"
TKF2[["percent.mt"]] <- PercentageFeatureSet(TKF2, pattern =  "^mt-")
TKF2[["percent.RPS"]] <- PercentageFeatureSet(TKF2, pattern = "^Rps")
TKF2[["percent.RPL"]] <- PercentageFeatureSet(TKF2, pattern = "^Rpl")
TKF2 <- subset(TKF2, subset = percent.mt <50)
TKF2 <- subset(TKF2, subset = percent.RPS <25)

#Foxd1;DLX HOMO
TKD8 <- Read10X_h5("HyDD2-HD2/scRNA/TKD8/raw_feature_bc_matrix.h5")
colnames(TKD8) = paste0("TKD8_", colnames(TKD8))
TKD8 <- CreateSeuratObject(counts = TKD8, project = "DLX",
                           min.cells = 5, min.features = 1000) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TKD8 <- subset(TKD8, subset = nCount_RNA > 2000)  
TKD8 <- RenameIdents(TKD8, "TKD8"= "Mut1")
TKD8 <- AddMetaData(TKD8, TKD8@active.ident, "Genotype")
Idents(TKD8) <- "Genotype"
TKD8[["percent.mt"]] <- PercentageFeatureSet(TKD8, pattern =  "^mt-")
TKD8[["percent.RPS"]] <- PercentageFeatureSet(TKD8, pattern = "^Rps")
TKD8[["percent.RPL"]] <- PercentageFeatureSet(TKD8, pattern =  "^Rpl")
TKD8 <- subset(TKD8, subset = percent.mt <50)
TKD8 <- subset(TKD8, subset = percent.RPS <25)

#Foxd1;DLX HOMO
TKF3 <- Read10X_h5("HyDD2-HD2/scRNA/TKF3/raw_feature_bc_matrix.h5")
colnames(TKF3) = paste0("TKF3_", colnames(TKF3))
TKF3 <- CreateSeuratObject(counts = TKF3, project = "DLX",
                           min.cells = 5, min.features = 1000) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TKF3 <- subset(TKF3, subset = nCount_RNA > 2000)  
TKF3 <- RenameIdents(TKF3, "TKF3"= "Mut2")
TKF3 <- AddMetaData(TKF3, TKF3@active.ident, "Genotype")
Idents(TKF3) <- "Genotype"
TKF3[["percent.mt"]] <- PercentageFeatureSet(TKF3, pattern = "^mt-")
TKF3[["percent.RPS"]] <- PercentageFeatureSet(TKF3, pattern = "^Rps")
TKF3[["percent.RPL"]] <- PercentageFeatureSet(TKF3, pattern = "^Rpl")
TKF3 <- subset(TKF3, subset = percent.mt <50)
TKF3 <- subset(TKF3, subset = percent.RPS <25)

#Foxd1;DLX HOMO
TKF4 <- Read10X_h5("HyDD2-HD2/scRNA/TKF4/raw_feature_bc_matrix.h5")
colnames(TKF4) = paste0("TKF4_", colnames(TKF4))
TKF4 <- CreateSeuratObject(counts = TKF4, project = "DLX",
                           min.cells = 5, min.features = 1000) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TKF4 <- subset(TKF4, subset = nCount_RNA > 2000)  
TKF4 <- RenameIdents(TKF4, "TKF4"= "Mut3")
TKF4 <- AddMetaData(TKF4, TKF4@active.ident, "Genotype")
Idents(TKF4) <- "Genotype"
TKF4[["percent.mt"]] <- PercentageFeatureSet(TKF4, pattern = "^mt-")
TKF4[["percent.RPS"]] <- PercentageFeatureSet(TKF4, pattern = "^Rps")
TKF4[["percent.RPL"]] <- PercentageFeatureSet(TKF4, pattern = "^Rpl")
TKF4 <- subset(TKF4, subset = percent.mt <50)
TKF4 <- subset(TKF4, subset = percent.RPS <25)

#checking values
mean(TKD7$nCount_RNA)
mean(TKD7$nFeature_RNA)
mean(TKF2$nCount_RNA)
mean(TKF2$nFeature_RNA)
mean(TKD8$nCount_RNA)
mean(TKD8$nFeature_RNA)
mean(TKF3$nCount_RNA)
mean(TKF3$nFeature_RNA)
mean(TKF4$nCount_RNA)
mean(TKF4$nFeature_RNA)

#merge all and clean up before merging to HyDD2 E12
DLXCKO <- merge(x = TKD7, y = list(TKF2, TKD8, TKF3, TKF4))
save(DLXCKO , file = "Data/scRNA/Robj/DLXCKO_P21_ZITRN.Robj")

rm(TKD7)
rm(TKF2)
rm(TKF3)
rm(TKF4)
rm(TKD8)

DLXCKO[["percent.sex"]] <- PercentageFeatureSet(DLXCKO,
                                                features = c("Xist","Malat1","Tsix"))
VlnPlot(DLXCKO, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                             "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)                 

#process
DLXCKO <- SCTransform(DLXCKO,
                      vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
DLXCKO<- RunPCA(DLXCKO, npcs = 50, ndims.print = NA, verbose = F)
DLXCKO<- RunHarmony(DLXCKO, group.by.vars = c("orig.ident"),
                    assay.use = "SCT", plot.convergence = T)
DLXCKO<- RunUMAP(DLXCKO, reduction = "harmony", dims = 1:20,
                 n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(DLXCKO, reduction = "umap", label = F,
        pt.size = 0.1, group.by = "Genotype")

Idents(DLXCKO) <- "Genotype"
DLXCKO <- RenameIdents(DLXCKO, "Ctrl1" = "Ctrl", "Ctrl2" = "Ctrl",
                       "Mut3" = "Mut",  "Mut1" = "Mut",  "Mut2" = "Mut")
DLXCKO <- AddMetaData(DLXCKO, DLXCKO@active.ident, "Genotype_Final")

FeaturePlot(DLXCKO, c("Aldh1l1","Mobp","Slc17a7","Slc32a1"), order = T)
VlnPlot(DLXCKO, c("Slc32a1","Gal"), group.by = "Genotype_Final")
VlnPlot(DLXCKO, c("Slc32a1","Gfap"), group.by = "Genotype_Final")

####Bring Atlas####
DLXCKO <- FindNeighbors(DLXCKO, dims = 1:20, reduction = "harmony")
DLXCKO <- FindClusters(DLXCKO, resolution = 0.8) #SCT_snn_res.0.8
DimPlot(DLXCKO, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

#GFAP
#MOBP higher

DLXCKO <- PrepSCTFindMarkers(DLXCKO)
markers <- FindAllMarkers(DLXCKO, test.use = "wilcox",
                          logfc.threshold = 1,
                          min.pct = 0.5, verbose = T)
write.csv(markers, file = "Data//DEG-DLXCKO.csv")

save(DLXCKO , file = "Data/scRNA/Robj/DLXCKO_P21_ZITRN.Robj") #leave non-neuronal cells

#0 = 
#1 = 
#2 = Glia
#3 = 
#4 = Glia 
#5 = Junk
#6 = Glia 
#7 = 
#8 = Glia  
#9 = Glia   
#10 = Glia  
#11 = 
#12 =  Glia 
#13 =  Glia 
#14 =  Glia 
#15 =  Glia 
#16 =  Glia 
#17 = 
#18 =  Glia 
#19 = 
#20 = 
#21 = 
#22 = 
#23 =  Glia 
#24 =  Glia 
#25 =  Glia 
#26 =  Glia 
#27 =
#28 =  Glia 
#29 = 
#30 =
#Remove non-neuronal cells 
DLXCKO <- subset(DLXCKO, idents = c("0", "1", "3", "7", 
                                    "11",  "17", "19", "20", 
                                    "21", "22",  "27", "29", "30"))
DLXCKO <- SCTransform(DLXCKO,
                      vars.to.regress = c("nCount_RNA","nFeature_RNA"))
DLXCKO <- RunPCA(DLXCKO, npcs = 50, ndims.print = NA, verbose = F)
DLXCKO<- RunHarmony(DLXCKO, group.by.vars = c("orig.ident"),
                    assay.use = "SCT", plot.convergence = T)
DLXCKO <- RunUMAP(DLXCKO, reduction = "harmony", dims = 1:20,
                  n.neighbours = 10L, min.dist = 0.01, spread = 3)
DimPlot(DLXCKO, reduction = "umap", label = T, pt.size = 0.1) + NoLegend() + NoAxes()

DLXCKO <- FindNeighbors(DLXCKO, dims = 1:20, reduction = "harmony")
DLXCKO <- FindClusters(DLXCKO, resolution = 0.8) #SCT_snn_res.0.8
DimPlot(DLXCKO, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

markers <- FindAllMarkers(DLXCKO, test.use = "wilcox",
                          logfc.threshold = 1,
                          min.pct = 0.5, verbose = T)
write.csv(markers, file = "Data//DEG-DLXCKO.csv")

FeaturePlot(DLXCKO, c("Slc17a7","Slc32a1"), order = T)
#second cleanup
#0 = 
#1 = 
#2 = 
#3 = 
#4 = 
#5 = 
#6 = 
#7 = 
#8 = 
#9 = 
#10 = 
#11 = 
#12 = Junk
#13 = 
#14 = Junk
#15 = Junk 
#16 = Junk  
#17 = 
#18 = Junk
#19 = 
#20 = 
#21 = 
#22 = Junk
#23 = Junk 
#24 = Junk


DLXCKO <- subset(DLXCKO, idents = c("0", "1", "2", "3", "4", "5", 
                                    "6", "7", "8", "9", "10", 
                                    "11", "13", "17",  "19", "20", 
                                    "21"))
DimPlot(DLXCKO, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
DLXCKO <- SCTransform(DLXCKO,
                      vars.to.regress = c("nCount_RNA","nFeature_RNA"))
DLXCKO <- RunPCA(DLXCKO, npcs = 50, ndims.print = NA, verbose = F)
DLXCKO<- RunHarmony(DLXCKO, group.by.vars = c("orig.ident"),
                    assay.use = "SCT", plot.convergence = T)
DLXCKO <- RunUMAP(DLXCKO, reduction = "harmony", dims = 1:20,
                  n.neighbours = 10L, min.dist = 0.01, spread = 3)
DimPlot(DLXCKO, reduction = "umap", label = T, pt.size = 0.1) + NoLegend() + NoAxes()

DLXCKO <- FindNeighbors(DLXCKO, dims = 1:20, reduction = "harmony")
DLXCKO <- FindClusters(DLXCKO, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(DLXCKO, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

markers <- FindAllMarkers(DLXCKO, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.2, verbose = T)
write.csv(markers, file = "Data//DEG-DLXCKO.csv")
save(DLXCKO , file = "Data/scRNA/Robj/DLXCKO_P21_ZITRN_Neuron.Robj")


#second cleanup
VlnPlot(DLXCKO, c("Slc32a1","Slc17a6","Slc17a7")) + NoLegend()
FeaturePlot(DLXCKO, c("Slc32a1","Slc17a7","Slc17a6"))

FeaturePlot(DLXCKO, c("Gal","Slc17a7"), blend = T, order = T, blend.threshold = 1)

cells <- WhichCells(DLXCKO, expression = Gal > 1 & Slc17a7 > 1)
cells <- WhichCells(DLXCKO, expression = Gal > 1 & Slc17a6 > 1)
DimPlot(DLXCKO, cells.highlight = cells)

FeaturePlot(DLXCKO, c("Gal","Slc17a6"), blend = T, order = T)

#0 = UNKNOWN, Glutamatergic (Slc17a6)
#1 = Lef1, Glutamatergic (Slc17a6), Glutamatergic (Slc17a7)
#2 = Calb1, Glutamatergic (Slc17a6)
#3 = Sst, GABAergic, Glutamatergic (Slc17a6)
#4 = Meis2, GABAergic
#5 = Lhx1os, GABAergic, Glutamatergic (Slc17a6)
#6 = Nfib, GABAergic, Glutamatergic (Slc17a6)
#7 = Dlx6os1, GABAergic
#8 = Pax6 / Meis1, GABAergic
#9 = Oxt / Esr2 / Sim1, Glutamatergic (Slc17a6)
#10 = Hcrt / LHx9, Glutamatergic (Slc17a6)
#11 = Crh / Avp / Trh, Glutamatergic (Slc17a6)
#12 = Foxg1 / Foxo1
#13 = Meis1 / Tbr1 / Lhx1, Glutamatergic (Slc17a6)
#14 = Pmch / Otx1 , Glutamatergic (Slc17a6)
#15 = Remove, Glutamatergic (Slc17a6), Glutamatergic (Slc17a7)
#16 = Remove, Glutamatergic (Slc17a6)

DLXCKO <- subset(DLXCKO, idents = c("0", "1", "2", "3", 
                                    "4", "5", "6", "7", "8",
                                    "9", "10", "11", "12",
                                    "13", "14"))
DLXCKO <- RenameIdents(DLXCKO, "0" = "Glutamatergic", "1" = "Lef1+",
                       "2" = "Calb1+", "3" = "Sst+", 
                       "4" = "Meis2+", "5" = "Lhx1os+",
                       "6" = "Nfib+", "7" = "Dlx6os1+",
                       "8" = "Pax6+", "9" = "Oxt+",
                       "10" = "Hcrt+", "11" = "Crh+",
                       "12" = "Foxg1+", "13" = "Tbr1+",
                       "14" = "Pmch+")
DLXCKO <- AddMetaData(DLXCKO, DLXCKO@active.ident, "Clusters")
DimPlot(DLXCKO, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

DLXCKO <- SCTransform(DLXCKO,
                      vars.to.regress = c("nCount_RNA","nFeature_RNA"))
DLXCKO <- RunPCA(DLXCKO, npcs = 50, ndims.print = NA, verbose = F)
DLXCKO<- RunHarmony(DLXCKO, group.by.vars = c("orig.ident"),
                    assay.use = "SCT", plot.convergence = T)
DLXCKO <- RunUMAP(DLXCKO, reduction = "harmony", dims = 1:20,
                  n.neighbours = 10L, min.dist = 0.01, spread = 3)
DimPlot(DLXCKO, reduction = "umap", label = T, pt.size = 0.1) + NoLegend() + NoAxes()
save(DLXCKO , file = "Data/scRNA/Robj/DLXCKO_P21_ZITRN_Neuron.Robj")


#match numbers
# Identify columns
Control_cols <- grep("^(TKD7_|TKF2)", colnames(DLXCKO), value = TRUE)
Mutant_cols <- grep("^(TKF3_|TKF4_|TKD8_)", colnames(DLXCKO), value = TRUE)

# Sample TK43_ columns
if(length(Mutant_cols) >= 17880) {
  sampled_Mutant_cols <- sample(Mutant_cols, size = 17880, replace = FALSE)
} else {
  sampled_Mutant_cols<- Mutant_cols
  message("Not enough columns to sample, using all available.")
}

# Combine with all columns
final_columns <- c(sampled_Mutant_cols, Control_cols)
# Subset the data frame
DLXCKO <- subset(DLXCKO, cells = final_columns)
DimPlot(DLXCKO, reduction = "umap", label = T,
        pt.size = 0.1, split.by = "Genotype_Final") + NoLegend() + NoAxes()
save(DLXCKO , file = "Data/scRNA/Robj/DLXCKO_P21_ZITRN_Neuron_Final.Robj")

####continue#####
markers <- FindAllMarkers(DLXCKO, test.use = "wilcox",
                          logfc.threshold = 0.2,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "Data/DEG-DLXCKO.csv")

Idents(DLXCKO) <- "Genotype_Final"
markers <- FindAllMarkers(DLXCKO, test.use = "wilcox",
                          logfc.threshold = 0.2,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "Data//DEG-DLXCKO_Genotype.csv")

####Density####
library(viridis)
library(viridisLite)
library(ggplot2)
Idents(DLXCKO) <- "Genotype_Final"
my_levels <- c("Ctrl","Mut")
Grid <- DLXCKO@reductions$umap@cell.embeddings
Grid <- data.frame(Grid)
Grid$orig.ident <- DLXCKO@meta.data$Genotype_Final
Grid$orig.ident = factor(Grid$orig.ident, levels = my_levels)
Grid_ds <- Grid[,1:2]
ggplot(Grid, aes(x = umap_1, y = umap_2)) +
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

####identify cluster####
Idents(Mutant) <- "predicted.id"
Cells1 <-WhichCells(Mutant, idents = "AntID_ID_TT")
Cells1 <-WhichCells(Mutant, idents = "SMN")
Cells1 <-WhichCells(Mutant, idents = "MMN")
DimPlot(Mutant, reduction = "umap", label = F, 
        pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
  ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf"))

####Percentage####
library(dplyr)
library(reshape2)
library(Rmisc)
meta.data <- DLXCKO@meta.data
Idents(DLXCKO) <- "Clusters"
graph <- 100*prop.table(table(Idents(DLXCKO), DLXCKO@meta.data$Genotype_Final), margin = 1) #margin 1= row, 2 = column
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




####CELL SELECT####
plot <- DimPlot(object = Mutant)
# Follow instructions in the terminal to select points
cells.located <- CellSelector(plot = plot)
cells.located
Idents(Mutant, cells = cells.located) <- "SelectedCells"
# Automatically set the identity class of selected cells and return a new Seurat object
Mutant <- CellSelector(plot = plot, object = Mutant, ident = 'SelectedCells')
table(Mutant@active.ident, Mutant@meta.data$Genotype)
DimPlot(Mutant)
test <- FindMarkers(Mutant, ident.1 = "SelectedCells",
                    logfc.threshold = 0.2,
                    min.pct = 0.05, verbose = T)

markers <- FindMarkers(Mutant,
                       test.use = "LR", ident.1 = "SelectedCells",
                       latent.vars = c("nCount_RNA","nFeature_RNA"),
                       logfc.threshold = 0.1,
                       min.pct = 0.1, verbose = T)

#WhichCells(Mutant, cells = cells.located)
#sub_obj <- subset(Mutant, cells = sub_cells)



####Junk####
#organize the order
categoryarray = c()
Mutant@active.ident = factor(Mutant@active.ident,
                             levels = categoryarray) 
Mutant@meta.data$Mutant_Broad = factor(Mutant@meta.data$Mutant_Broad,
                                       levels = categoryarray) 
#plot
DimPlot(Mutant) + NoAxes()





