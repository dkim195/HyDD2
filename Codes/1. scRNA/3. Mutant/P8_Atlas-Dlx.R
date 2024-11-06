####Bring Atlas####
library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("E:/")
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

#atlas
load(file = "Robj/P8_Neuron.Robj")
P8@meta.data$Cluster_Pass1a <- NULL
P8@meta.data$SCT_snn_res.2 <- NULL
P8@meta.data$seurat_clusters <- NULL
Idents(P8) <- "Age"
P8 <- RenameIdents(P8,"P8_Rep1" = "Atlas","P8_Rep2" = "Atlas","TKC5" = "Atlas","TKC6" = "Atlas")
P8 <- AddMetaData(P8, P8@active.ident, "Reference")
Idents(P8) <- "Cluster_Pass3"

#Lhx2
load(file = "Robj/P8_Dlx_Neuron.Robj")
P8_Dlx$SCT_snn_res.0.8 <- NULL
P8_Dlx$old.ident <- NULL
P8_Dlx$seurat_clusters <- NULL
Idents(P8_Dlx) <- "Age"
P8_Dlx <- RenameIdents(P8_Dlx, "P8_Dlx_Rep1" = "Mutant", "P8_Dlx_Rep2" = "Mutant" , "P8_Dlx_Rep3" = "Mutant")
P8_Dlx <- AddMetaData(P8_Dlx, P8_Dlx@active.ident, "Reference")
P8_Dlx <- AddMetaData(P8_Dlx, P8_Dlx@active.ident, "Cluster_Pass3")
Idents(P8_Dlx) <- "Cluster_Pass3"

####transfer method####
options(future.globals.maxSize = 8000 * 1024^2)

Mutant.anchors <- FindTransferAnchors(reference = P8, query = P8_Dlx,
                                      normalization.method = "SCT",
                                      reduction = "pcaproject")
predictions <- TransferData(anchorset = Mutant.anchors,
                            refdata = P8$Cluster_Pass3, dims = 1:30)
P8_Dlx <- AddMetaData(P8_Dlx, metadata = predictions)

P8_Dlx$prediction.match <- P8_Dlx$predicted.id == P8_Dlx$orig.ident
table(P8_Dlx$prediction.match) ##0 as PC0 is similar as PR0
DimPlot(P8_Dlx, reduction = "umap", group.by = "predicted.id", label=T)
P8_Dlx@meta.data$Cluster_Pass3 <- P8_Dlx@meta.data$predicted.id

####final####
#follows after transfer
Idents(P8_Dlx) <- "Cluster_Pass3"

Idents(P8) <- "Cluster_Pass3"
Mutant <- merge(x = P8, y = list(P8_Dlx))
Mutant[["percent.sex"]] <- PercentageFeatureSet(Mutant,
                                                features = c("Xist","Malat1","Tsix"))
VlnPlot(Mutant, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                             "percent.RPS","percent.RPL","percent.sex"), pt.size = 0, group.by = "Reference")  
save(Mutant, file = "Robj/P8_Dlx_Final.Robj")

####Final_Processing####
library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
setwd("E:/")
load(file = "Robj/P8_Dlx_Final.Robj")
Mutant <- SCTransform(Mutant, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
Mutant<- RunPCA(Mutant, npcs = 50, ndims.print = NA, verbose = F)
Mutant<- RunHarmony(Mutant, group.by.vars = c("orig.ident"), assay.use = "SCT", plot.convergence = T)
Mutant<- RunUMAP(Mutant, reduction = "harmony", dims = 1:22,
                 n.neighbours = 10L, min.dist = 0.01, spread = 1)
Idents(Mutant) <- "Reference"
Mutant <- RenameIdents(Mutant, "Atlas" = "Control")
Mutant <- AddMetaData(Mutant, Mutant@active.ident, "Genotype")

#DefaultAssay(Mutant) <- "RNA"
FeaturePlot(Mutant, "Dlx1", split.by = "Genotype", order = T )

DimPlot(Mutant, reduction = "umap", label = F,
        pt.size = 0.1, split.by = "Genotype") + NoLegend() + NoAxes()

DimPlot(Mutant, reduction = "umap", label = T,
        group.by = "Cluster_Pass3",split.by = "Genotype") + NoLegend() + NoAxes()

DimPlot(Mutant, reduction = "umap", label = T,
        group.by = "Cluster_Pass3") + NoLegend() + NoAxes()

VlnPlot(Mutant, c("Nr4a2","Foxa1","Irx5","Foxb1","Lhx6","Sp9","Pitx2"), group.by = "Genotype")

save(Mutant, file = "Robj/P8_Dlx_Final.Robj")



FeaturePlot(Mutant, "Sp8", split.by = "Genotype", order = T )
FeaturePlot(Mutant, "Irx5", split.by = "Genotype", order = T )
FeaturePlot(Mutant, "Pitx2", split.by = "Genotype", order = T )
FeaturePlot(Mutant, "Meis2", split.by = "Genotype", order = T )
FeaturePlot(Mutant, "Isl1", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Rax", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Sim1", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Vax1", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Six3", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Lhx5", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Lhx1", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Arx", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Rorb", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Sst", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Avp", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Foxg1", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Cck", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Trh", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Tcf7l2", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Hmx2", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Gal", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Npy", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Calb1", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Cck", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Meis2", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Slc32a1", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Gad1", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Gad2", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Pax6", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Sp8", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Sp9", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Dlx6os1", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Dlx2", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Dlx1", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Lhx6", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Vipr2", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Cdh23", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Oxt", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Slc6a3", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Th", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Ddc", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Crh", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Onecut2", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Onecut3", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Neurod2", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Tbr1", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Eomes", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Satb2", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Emx2", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Zbtb18", split.by = "Genotype") + NoLegend()
FeaturePlot(Mutant, "Kdm5d", split.by = "Genotype") + NoLegend()

VlnPlot(Mutant, c("Avp","Gal","Calb1","Oxt",
                  "Meis2","Gad1","Slc32a1","Slc17a6"), group.by = "Genotype", pt.size = 0.1) + NoLegend()

VlnPlot(Mutant, c("Cck","Lhx6","Pvalb","Sst",
                  "Tac1","Tac2","Npy","Nts"), group.by = "Genotype", pt.size = 0.1) + NoLegend()

VlnPlot(Mutant, c("Pnoc","Dlx1","Dlx2","Sp8",
                  "Nfia","Hdc","Gal","Th"), group.by = "Genotype", pt.size = 0.1) + NoLegend()

VlnPlot(Mutant, c("Agrp","Pomc","Hcrt","Pmch",
                  "Kiss1","Ghrh","Npw","Grp"), group.by = "Genotype", pt.size = 0.1) + NoLegend()

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
Idents(Mutant) <- "Cluster_Pass3"
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

#WhichCells(Mutant, cells = cells.located)
#sub_obj <- subset(Mutant, cells = sub_cells)

