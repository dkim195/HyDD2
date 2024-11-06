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

####Nkx2.2 Cre-knock in####
#E12
#10x v2
setwd("D:/")

#TK33
TK33 <- Read10X(data.dir = "scRNA/TK33/raw_gene_bc_matrices/mm10/")
colnames(TK33) = paste0("TK33_", colnames(TK33))
TK33 <- CreateSeuratObject(counts = TK33, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK33 <- subset(TK33, subset = nCount_RNA > 1000)  
TK33 <- RenameIdents(TK33, "TK33"= "Nkx22_Mut1")
TK33 <- AddMetaData(TK33, TK33@active.ident, "Age")
Idents(TK33) <- "Age"
TK33[["percent.mt"]] <- PercentageFeatureSet(TK33, pattern = "^mt-")
TK33[["percent.RPS"]] <- PercentageFeatureSet(TK33, pattern = "^Rps")
TK33[["percent.RPL"]] <- PercentageFeatureSet(TK33, pattern = "^Rpl")
TK33 <- subset(TK33, subset = percent.mt <50)
TK33 <- subset(TK33, subset = percent.RPS <25)

#TK34
TK34 <- Read10X(data.dir = "scRNA/TK34/raw_gene_bc_matrices/mm10/")
colnames(TK34) = paste0("TK34_", colnames(TK34))
TK34 <- CreateSeuratObject(counts = TK34, project = "E12",
                           min.cells = 5, min.features = 500) #1000 if depth is higher AND THEN CHANGE UMI to 2000
TK34 <- subset(TK34, subset = nCount_RNA > 1000)  
TK34 <- RenameIdents(TK34, "TK34"= "Nkx22_Mut2")
TK34 <- AddMetaData(TK34, TK34@active.ident, "Age")
Idents(TK34) <- "Age"
TK34[["percent.mt"]] <- PercentageFeatureSet(TK34, pattern = "^mt-")
TK34[["percent.RPS"]] <- PercentageFeatureSet(TK34, pattern = "^Rps")
TK34[["percent.RPL"]] <- PercentageFeatureSet(TK34, pattern = "^Rpl")
TK34 <- subset(TK34, subset = percent.mt <50)
TK34 <- subset(TK34, subset = percent.RPS <25)

#merge all and clean up before merging to HyDD2 E12
Nkx22 <- merge(x = TK33, y = list(TK34))
Nkx22[["percent.sex"]] <- PercentageFeatureSet(Nkx22,
                                                   features = c("Xist","Malat1","Tsix"))
VlnPlot(Nkx22, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                                "percent.RPS","percent.RPL","percent.sex"), pt.size = 0)                 

#process
Nkx22 <- SCTransform(Nkx22, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
Nkx22<- RunPCA(Nkx22, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(Nkx22, ndims = 50, reduction = "pca" )
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene)
m.s.genes <- homologene::human2mouse(s.genes,
                                     db = homologeneData2) 
m.g2m.genes <- homologene::human2mouse(g2m.genes,
                                       db = homologeneData2) 
Nkx22 <- CellCycleScoring(Nkx22, s.features = m.s.genes$mouseGene,
                              g2m.features = m.g2m.genes$mouseGene, set.ident = TRUE)
Nkx22<- RunHarmony(Nkx22, group.by.vars = c("orig.ident","Phase"), assay.use = "SCT", plot.convergence = T)
Nkx22<- RunUMAP(Nkx22, reduction = "harmony", dims = 1:20,
                    n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(Nkx22, reduction = "umap", label = F,pt.size = 0.1, group.by = "Phase")
DimPlot(Nkx22, reduction = "umap", label = F,pt.size = 0.1, group.by = "Age")
Idents(Nkx22) <- "Age"
DimPlot(Nkx22, reduction = "umap", label = F,pt.size = 0.1, cols = , 
        split.by = "Age", ncol = 2) + NoLegend()
DimPlot(Nkx22, reduction = "umap", label = F,pt.size = 0.1, cols =, group.by = "Age", 
        split.by = "Phase", ncol = 2) + NoLegend()

FeaturePlot(Nkx22, c("Sp8","Pitx2","Meis2","Sp9", "Irx5"))
save(Nkx22, file = "Robj/Nkx22.Robj")


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

load(file = "Robj/Nkx22.Robj")
Idents(Nkx22) <- "Age"
Nkx22 <- RenameIdents(Nkx22, "Nkx22_Mut1" = "Mutant", "Nkx22_Mut2" = "Mutant")
Nkx22 <- AddMetaData(Nkx22, Nkx22@active.ident, "Reference")
Nkx22 <- AddMetaData(Nkx22, Nkx22@active.ident, "Cluster_Pass2")
Idents(Nkx22) <- "Cluster_Pass2"


options(future.globals.maxSize = 8000 * 1024^2)
E12_Atlas <- E12_Atlas[, sample(colnames(E12_Atlas), size = 72930, replace=F)] #70000 seems to work

Mutant.anchors <- FindTransferAnchors(reference = E12_Atlas, query = Nkx22,
                                      normalization.method = "SCT",
                                      reduction = "pcaproject")
predictions <- TransferData(anchorset = Mutant.anchors,
                            refdata = E12_Atlas$Cluster_Pass2, dims = 1:30)
Nkx22 <- AddMetaData(Nkx22, metadata = predictions)
DimPlot(Nkx22, group.by = "predicted.id")

Nkx22$prediction.match <- Nkx22$predicted.id == Nkx22$orig.ident
table(Nkx22$prediction.match) ##0 as PC0 is similar as PR0
DimPlot(Nkx22, reduction = "umap", group.by = "predicted.id")
Nkx22@meta.data$Cluster_Pass2 <- Nkx22@meta.data$predicted.id

####final####
#follows after transfer
Idents(Nkx22) <- "Cluster_Pass2"
E12_Atlas <- subset(E12_Atlas, subset = nCount_RNA < 10000) #reduce to match quality for v2
E12_Atlas <- subset(E12_Atlas, subset = nFeature_RNA < 2000) #reduce to match quality for v2
E12_Atlas <- E12_Atlas[, sample(colnames(E12_Atlas), size =7959, replace=F)] #size = cell no similar to Nkx22
Idents(E12_Atlas) <- "Cluster_Pass2"
Mutant <- merge(x = E12_Atlas, y = list(Nkx22))
Mutant[["percent.sex"]] <- PercentageFeatureSet(Mutant,
                                                features = c("Xist","Malat1","Tsix"))
VlnPlot(Mutant, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                             "percent.RPS","percent.RPL","percent.sex"), pt.size = 0, group.by = "Reference")                 
Mutant <- SCTransform(Mutant, vars.to.regress = c("nCount_RNA","nFeature_RNA")) #CHECK SEX AND MITOCHODNRIA
Mutant<- RunPCA(Mutant, npcs = 50, ndims.print = NA, verbose = F)
ElbowPlot(Mutant, ndims = 50, reduction = "pca" )
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

save(Mutant, file = "Robj/Nkx22_Final.Robj")


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





