library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)

#install.packages('https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz',
#                 repos=NULL,type="source", INSTALL_opts = "--no-lock")
#or remotes::install_version("spatstat", version = "1.64-1") space 
#remotes::install_version("Seurat", version = "3.2.3")
#install.packages("htmltools") double check 
#ggplot2 v3.2.1 -> ggplot v3.3.5 using install.packages("ggplot2")
#install.packages("sctransform") to 0.3.2 to update
#install.packages("uwot")
#the version i have been using is remotes::install_version("Seurat", version = "3.1.5")
#check https://github.com/satijalab/seurat/releases

####Foxd1 Isl1####
#E12
#10x v3
setwd("D:/")

#TK37 # low quality
TK37 <- Read10X_h5("scRNA/TK37/outs/raw_feature_bc_matrix.h5")
colnames(TK37) = paste0("TK37_", colnames(TK37))
TK37 <- CreateSeuratObject(counts = TK37, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK37 <- subset(TK37, subset = nCount_RNA > 1000)  
TK37 <- RenameIdents(TK37, "TK37"= "Isl1_Mut1")
TK37 <- AddMetaData(TK37, TK37@active.ident, "Age")
Idents(TK37) <- "Age"
TK37[["percent.mt"]] <- PercentageFeatureSet(TK37, pattern = "^mt-")
TK37[["percent.RPS"]] <- PercentageFeatureSet(TK37, pattern = "^Rps")
TK37[["percent.RPL"]] <- PercentageFeatureSet(TK37, pattern = "^Rpl")
TK37 <- subset(TK37, subset = percent.mt <50)
TK37 <- subset(TK37, subset = percent.RPS <25)

#TK39
TK39 <- Read10X_h5("scRNA/TK39/outs/raw_feature_bc_matrix.h5")
colnames(TK39) = paste0("TK39_", colnames(TK39))
TK39 <- CreateSeuratObject(counts = TK39, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK39 <- subset(TK39, subset = nCount_RNA > 1000)  
TK39 <- RenameIdents(TK39, "TK39"= "Isl1_Mut2")
TK39 <- AddMetaData(TK39, TK39@active.ident, "Age")
Idents(TK39) <- "Age"
TK39[["percent.mt"]] <- PercentageFeatureSet(TK39, pattern = "^mt-")
TK39[["percent.RPS"]] <- PercentageFeatureSet(TK39, pattern = "^Rps")
TK39[["percent.RPL"]] <- PercentageFeatureSet(TK39, pattern = "^Rpl")
TK39 <- subset(TK39, subset = percent.mt <50)
TK39 <- subset(TK39, subset = percent.RPS <25)

#TK41 - METHANOL FIXED, higher threshold is needed
TK41 <- Read10X_h5("scRNA/TK41/outs/raw_feature_bc_matrix.h5")
colnames(TK41) = paste0("TK41_", colnames(TK41))
TK41 <- CreateSeuratObject(counts = TK41, project = "E12",
                           min.cells = 5, min.features = 1000) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK41 <- subset(TK41, subset = nCount_RNA > 2000)  
TK41 <- RenameIdents(TK41, "TK41"= "Isl1_Mut3")
TK41 <- AddMetaData(TK41, TK41@active.ident, "Age")
Idents(TK41) <- "Age"
TK41[["percent.mt"]] <- PercentageFeatureSet(TK41, pattern = "^mt-")
TK41[["percent.RPS"]] <- PercentageFeatureSet(TK41, pattern = "^Rps")
TK41[["percent.RPL"]] <- PercentageFeatureSet(TK41, pattern = "^Rpl")
TK41 <- subset(TK41, subset = percent.mt <50)
TK41 <- subset(TK41, subset = percent.RPS <25)

mean(TK37$nCount_RNA)
mean(TK37$nFeature_RNA)
mean(TK39$nCount_RNA)
mean(TK39$nFeature_RNA)
mean(TK41$nCount_RNA)
mean(TK41$nFeature_RNA)

#merge all and clean up before merging to HyDD2 E12
Isl1 <- merge(x = TK37, y = list(TK39, TK41))
Isl1[["percent.sex"]] <- PercentageFeatureSet(Isl1,
                                               features = c("Xist","Malat1","Tsix"))
VlnPlot(Isl1, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                            "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)                 
#check cell no and quality and identify the parameter to adjust the Atlas
#process
Isl1 <- SCTransform(Isl1, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
Isl1<- RunPCA(Isl1, npcs = 50, ndims.print = NA, verbose = F)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
Isl1 <- CellCycleScoring(Isl1, s.features = m.s.genes$mouseGene,
                          g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
Isl1<- RunHarmony(Isl1, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
Isl1<- RunUMAP(Isl1, reduction = "harmony", dims = 1:20,
                n.neighbours = 10L, min.dist = 0.01, spread = 1)
save(Isl1, file = "Robj/Isl1.Robj")
DimPlot(Isl1, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(Isl1, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
Idents(Isl1) <- "Age"
DimPlot(Isl1, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(Isl1, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
        split.by = "Phase", ncol = 2) + NoLegend()

FeaturePlot(Isl1, c("Sp8","Pitx2","Meis2","Sp9", "Irx5"))

####Bring Atlas####
library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("D:/")
load(file = "Robj/E12_Atlas_Merge_V2.Robj")
E12_Atlas@meta.data$old.ident <- NULL
E12_Atlas@meta.data$Cluster_Pass1 <- NULL
E12_Atlas@meta.data$Cluster <- NULL
E12_Atlas@meta.data$SCT_snn_res.0.4 <- NULL
E12_Atlas@meta.data$seurat_clusters <- NULL
Idents(E12_Atlas) <- "Age"

E12_Atlas <- RenameIdents(E12_Atlas, "E12" = "Atlas", "E12_Rep1" = "Atlas", "E12_Rep2" = "Atlas",
                          "E12_Atlas1" = "Atlas", "E12_Atlas2" = "Atlas", "E12_Atlas3" = "Atlas",
                          "E12_Atlas4" = "Atlas", "E12_Atlas5" = "Atlas", "E12_Atlas6" = "Atlas",
                          "E12_Atlas7" = "Atlas", "E12_Atlas8" = "Atlas", "E12_Atlas9" = "Atlas" )
E12_Atlas <- AddMetaData(E12_Atlas, E12_Atlas@active.ident, "Reference")
Idents(E12_Atlas) <- "Cluster_Pass2"

load(file = "Robj/Isl1.Robj")
Idents(Isl1) <- "Age"
Isl1 <- RenameIdents(Isl1, "Isl1_Mut1" = "Mutant", "Isl1_Mut2" = "Mutant", "Isl1_Mut3" = "Mutant")
Isl1 <- AddMetaData(Isl1, Isl1@active.ident, "Reference")
Isl1 <- AddMetaData(Isl1, Isl1@active.ident, "Cluster_Pass2")
Idents(Isl1) <- "Cluster_Pass2"

options(future.globals.maxSize = 8000 * 1024^2)

E12_Atlas <- E12_Atlas[, sample(colnames(E12_Atlas), size = 72930, replace=F)] #70000 seems to work


Mutant.anchors <- FindTransferAnchors(reference = E12_Atlas, query = Isl1,
                                      normalization.method = "SCT",
                                      reduction = "pcaproject")
predictions <- TransferData(anchorset = Mutant.anchors,
                            refdata = E12_Atlas$Cluster_Pass2, dims = 1:30)
Isl1 <- AddMetaData(Isl1, metadata = predictions)
Isl1$prediction.match <- Isl1$predicted.id == Isl1$orig.ident
table(Isl1$prediction.match) ##0 as PC0 is similar as PR0
DimPlot(Isl1, reduction = "umap", group.by = "predicted.id", label=T)
Isl1@meta.data$Cluster_Pass2 <- Isl1@meta.data$predicted.id

####final####
#follows after transfer
Idents(Isl1) <- "Cluster_Pass2"
E12_Atlas <- subset(E12_Atlas, subset = nCount_RNA < 10000) #reduce to match quality for v2
E12_Atlas <- subset(E12_Atlas, subset = nFeature_RNA < 3000) #reduce to match quality for v2

E12_Atlas <- E12_Atlas[, sample(colnames(E12_Atlas), size =21393, replace=F)] #size = cell no similar to Isl1
Idents(E12_Atlas) <- "Cluster_Pass2"
Mutant <- merge(x = E12_Atlas, y = list(Isl1))
Mutant[["percent.sex"]] <- PercentageFeatureSet(Mutant,
                                                features = c("Xist","Malat1","Tsix"))
VlnPlot(Mutant, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                             "percent.RPS","percent.RPL","percent.sex"), pt.size = 0, group.by = "Reference")                 
Mutant <- SCTransform(Mutant, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
Mutant<- RunPCA(Mutant, npcs = 50, ndims.print = NA, verbose = F)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
Mutant <- CellCycleScoring(Mutant, s.features = m.s.genes$mouseGene,
                           g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
Mutant<- RunHarmony(Mutant, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
Mutant<- RunUMAP(Mutant, reduction = "harmony", dims = 1:22,
                 n.neighbours = 10L, min.dist = 0.01, spread = 1)
Idents(Mutant) <- "Reference"
Mutant <- RenameIdents(Mutant, "Atlas" = "Control")
Mutant <- AddMetaData(Mutant, Mutant@active.ident, "Genotype")

#DefaultAssay(Mutant) <- "RNA"
FeaturePlot(Mutant, "Sp8", split.by = "Genotype", order = T )
FeaturePlot(Mutant, "Irx5", split.by = "Genotype", order = T )
FeaturePlot(Mutant, "Pitx2", split.by = "Genotype", order = T )
FeaturePlot(Mutant, "Meis2", split.by = "Genotype", order = T )

DimPlot(Mutant, reduction = "umap", label = F,
        pt.size = 0.1, split.by = "Genotype") + NoLegend() + NoAxes()

DimPlot(Mutant, reduction = "umap", label = T,
        group.by = "Cluster_Pass2",split.by = "Genotype") + NoLegend() + NoAxes()

DimPlot(Mutant, reduction = "umap", label = T,
        group.by = "Cluster_Pass2") + NoLegend() + NoAxes()

VlnPlot(Mutant, c("Nr4a2","Foxa1","Irx5","Foxb1","Lhx6","Sp9","Pitx2","Chchd10"), group.by = "Genotype")
save(Mutant, file = "Robj/Isl1_Final.Robj")

####Density####
library(viridis)
library(viridisLite)
library(ggplot2)
Idents(Mutant) <- "Genotype"
my_levels <- c("Control","Mutant")
Grid <- Mutant@reductions$umap@cell.embeddings
Grid <- data.frame(Grid)
Grid$orig.ident <- Mutant@meta.data$Genotype
Grid$orig.ident = factor(Grid$orig.ident, levels = my_levels)
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
meta.data <- Mutant@meta.data
Idents(Mutant) <- "Cluster_Pass2"
graph <- 100*prop.table(table(Idents(Mutant), Mutant@meta.data$Genotype), margin = 1) #margin 1= row, 2 = column
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

####Junk####
#organize the order
categoryarray = c()
Mutant@active.ident = factor(Mutant@active.ident,
                             levels = categoryarray) 
Mutant@meta.data$Mutant_Broad = factor(Mutant@meta.data$Mutant_Broad,
                                       levels = categoryarray) 
#plot
DimPlot(Mutant) + NoAxes()





