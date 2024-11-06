library(cowplot)
library(dplyr)
library(Matrix)
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(future)
library(harmony)
set.seed(1234)
plan("multicore", workers = 10)
plan()
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram

setwd("/media/thomaskim/Data/")

####Loading####
#ATAC

#control
load(file = "scATAC/Files/TK43.Robj")
Idents(TK43) <- "orig.ident"
TK43 <- RenameIdents(TK43, "ATAC" = "E12")
TK43 <- AddMetaData(TK43, TK43@active.ident, "Age")
Idents(TK43) <- "orig.ident"
TK43 <- RenameIdents(TK43, "ATAC" = "TK43")
TK43 <- AddMetaData(TK43, TK43@active.ident, "Sample")
Idents(TK43) <- "orig.ident"
TK43 <- RenameIdents(TK43, "ATAC" = "Control")
TK43 <- AddMetaData(TK43, TK43@active.ident, "Genotype")
Idents(TK43) <- "Age"

load(file = "scATAC/Files/TK47.Robj")
Idents(TK47) <- "orig.ident"
TK47 <- RenameIdents(TK47, "ATAC" = "E12")
TK47 <- AddMetaData(TK47, TK47@active.ident, "Age")
Idents(TK47) <- "orig.ident"
TK47 <- RenameIdents(TK47, "ATAC" = "TK47")
TK47 <- AddMetaData(TK47, TK47@active.ident, "Sample")
Idents(TK47) <- "orig.ident"
TK47 <- RenameIdents(TK47, "ATAC" = "Control")
TK47 <- AddMetaData(TK47, TK47@active.ident, "Genotype")
Idents(TK47) <- "Age"

#mutant
load(file = "scATAC/Files/TK57.Robj")
Idents(TK57) <- "orig.ident"
TK57 <- RenameIdents(TK57, "ATAC" = "E12")
TK57 <- AddMetaData(TK57, TK57@active.ident, "Age")
Idents(TK57) <- "orig.ident"
TK57 <- RenameIdents(TK57, "ATAC" = "TK57")
TK57 <- AddMetaData(TK57, TK57@active.ident, "Sample")
Idents(TK57) <- "orig.ident"
TK57 <- RenameIdents(TK57, "ATAC" = "Mutant")
TK57 <- AddMetaData(TK57, TK57@active.ident, "Genotype")
Idents(TK57) <- "Age"

load(file = "scATAC/Files/TK58.Robj")
Idents(TK58) <- "orig.ident"
TK58 <- RenameIdents(TK58, "ATAC" = "E12")
TK58 <- AddMetaData(TK58, TK58@active.ident, "Age")
Idents(TK58) <- "orig.ident"
TK58 <- RenameIdents(TK58, "ATAC" = "TK58")
TK58 <- AddMetaData(TK58, TK58@active.ident, "Sample")
Idents(TK58) <- "orig.ident"
TK58 <- RenameIdents(TK58, "ATAC" = "Mutant")
TK58 <- AddMetaData(TK58, TK58@active.ident, "Genotype")
Idents(TK58) <- "Age"

####Fragment####
DefaultAssay(TK43) <- 'peaks'
DefaultAssay(TK47) <- 'peaks'
DefaultAssay(TK57) <- 'peaks'
DefaultAssay(TK58) <- 'peaks'

frags <- Fragments(TK43)  # get list of fragment objects
Fragments(TK43) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK43/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK43) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK47)  # get list of fragment objects
Fragments(TK47) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK47/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK47) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK57)  # get list of fragment objects
Fragments(TK57) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK57/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK57) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK58)  # get list of fragment objects
Fragments(TK58) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK58/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK58) <- frags  # assign update list of fragment objects back to the assay

####Intersection####
# Find the union of peaks between the two datasets
# Extract peak regions from both datasets
DefaultAssay(TK43) <- 'peaks'
DefaultAssay(TK47) <- 'peaks'
DefaultAssay(TK57) <- 'peaks'
DefaultAssay(TK58) <- 'peaks'

peaks1 <- GetAssayData(TK43, assay = "peaks", slot = "ranges")
peaks2 <- GetAssayData(TK47, assay = "peaks", slot = "ranges")
peaks3 <- GetAssayData(TK57, assay = "peaks", slot = "ranges")
peaks4 <- GetAssayData(TK58, assay = "peaks", slot = "ranges")

# Find the intersection of peaks
combined.peaks <- reduce(x = c(peaks1, peaks2, peaks3, peaks4))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# Subset both objects to include only common peaks
DefaultAssay(TK43) <- 'peaks'
DefaultAssay(TK47) <- 'peaks'
DefaultAssay(TK57) <- 'peaks'
DefaultAssay(TK58) <- 'peaks'

TK43.counts <- FeatureMatrix(
  fragments = Fragments(TK43),
  features = combined.peaks,
  cells = colnames(TK43))
TK47.counts <- FeatureMatrix(
  fragments = Fragments(TK47),
  features = combined.peaks,
  cells = colnames(TK47))
TK57.counts <- FeatureMatrix(
  fragments = Fragments(TK57),
  features = combined.peaks,
  cells = colnames(TK57))
TK58.counts <- FeatureMatrix(
  fragments = Fragments(TK58),
  features = combined.peaks,
  cells = colnames(TK58))

TK43[['peakunion']] <- CreateAssayObject(counts = TK43.counts)
TK47[['peakunion']] <- CreateAssayObject(counts = TK47.counts)
TK57[['peakunion']] <- CreateAssayObject(counts = TK57.counts)
TK58[['peakunion']] <- CreateAssayObject(counts = TK58.counts)

DefaultAssay(TK43) <- 'peakunion'
DefaultAssay(TK47) <- 'peakunion'
DefaultAssay(TK57) <- 'peakunion'
DefaultAssay(TK58) <- 'peakunion'

rm(TK43.counts)
rm(TK47.counts)
rm(TK57.counts)
rm(TK58.counts)

####Merging####
# Merge the datasets
Mutant <- merge(TK43, y = c(TK47, TK57, TK58),
                add.cell.ids = c("TK43", "TK47", "TK57", "TK58"))

rm(TK43)
rm(TK47)
rm(TK57)
rm(TK58)

save(Mutant, file = "/media/thomaskim/Data/scATAC/Files/Dlx_Mutant_adjusted.Robj")

####ATAC Process####
# Run TF-IDF normalization
Mutant <- RunTFIDF(Mutant)

# Find variable features
Mutant <- FindTopFeatures(Mutant, min.cutoff = "q0")

# Run dimensional reduction (LSI)
Mutant <- RunSVD(
  object = Mutant,
  assay = 'peakunion', #peaks
  reduction.key = 'LSI_',
  reduction.name = 'lsi')

Mutant <- RunHarmony(object = Mutant,
                     assay.use = 'peaks',
                     group.by.vars = 'Sample',
                     reduction.use = 'lsi', project.dim = F)

# Run UMAP
#DepthCor(Mutant)
Mutant <- RunUMAP(Mutant, reduction = 'harmony', dims = 2:30)

# Find neighbors and clusters
Mutant <- FindNeighbors(Mutant, reduction = 'harmony', dims = 2:30)
#Mutant <- FindClusters(Mutant, verbose = FALSE)

# Plot UMAP
DimPlot(Mutant, label = , group.by = "Sample")
DimPlot(Mutant, label = , group.by = "Genotype")

####Transfer Label####
#RNA
load(file = "scRNA/Robj/E12_Atlas_Merge_V2.Robj")
E12_Atlas <- UpdateSeuratObject(E12_Atlas)
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

E12_Atlas <- RenameIdents(E12_Atlas, "Prethalamus (Sp8)" = "Prethalamus",
                          "Prethalamus (Sp9)" = "Prethalamus", "Prethalamus (Sst)" = "Prethalamus",
                          "PVH_SON (Cartpt)" = "PVH_SON", "PVH_SON (Otp)" = "PVH_SON")
E12_Atlas <- AddMetaData(E12_Atlas, E12_Atlas@active.ident, "Cluster_Pass2")
E12_Atlas <- E12_Atlas[, sample(colnames(E12_Atlas), size = 72930, replace=F)] #70000 seems to work
DefaultAssay(Mutant) <- 'RNA'
DefaultAssay(E12_Atlas) <- 'RNA'

E12_Atlas <- FindVariableFeatures(object = E12_Atlas,  nfeatures = 1000)

transfer.anchors <- FindTransferAnchors(
  reference = E12_Atlas, query = Mutant,
  normalization.method = "LogNormalize",
  reduction = "pcaproject", dims = 1:30)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = E12_Atlas$Cluster_Pass2,
  weight.reduction = Mutant[['lsi']],
  dims = 2:30)

Mutant <- AddMetaData(object = Mutant, metadata = predicted.labels)

DimPlot(object = Mutant, label = TRUE, group.by = "predicted.id") + NoLegend()
Idents(Mutant) <- "predicted.id"

FeaturePlot(
  object = Mutant,
  features = c("Sp8","Nr4a2","Pitx2","Hmx2","Hmx3"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3)

VlnPlot(Mutant, c("Gadd45g","Ascl1","Neurog2","Meis2",
                  "Otp","Hmx2","Nr4a2","Pax6","Mutant_adjusted","Foxb1"), same.y.lims = T,
        pt.size = 0, stack = T, flip = T)


library(dplyr)
library(reshape2)
library(Rmisc)
meta.data <- Mutant@meta.data
Idents(Mutant) <- "predicted.id"
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

rm(E12_Atlas)
save(Mutant, file = "/media/thomaskim/Data/scATAC/Files/Dlx_Mutant_adjusted.Robj")

Idents(Mutant) <- "Genotype"
DefaultAssay(Mutant) <- "peakunion"
da_peaksx <- FindAllMarkers(
  object = Mutant,
  test.use = 'wilcox')

VlnPlot(object = Mutant,
  features = rownames(da_peaks)[1],
  pt.size = 0.1)

open_l25 <- rownames(da_peaksx[da_peaksx$avg_log2FC > 3, ])
open_l457 <- rownames(da_peaksx[da_peaksx$avg_log2FC < 3, ])

DefaultAssay(Mutant) <- "peaks"
closest_l25 <- ClosestFeature(Mutant, regions = open_l25)
closest_l457 <- ClosestFeature(Mutant,regions =  open_l457)
#increase in Atoh1 and Atoh7

####Match cell number####
Mutant$Cluster_Pass2 <- Mutant$predicted.id
Idents(Mutant) <- "Genotype"

# Assuming 'Mutant' is your Seurat object
# Assume you have a metadata column identifying each cell as 'control' or 'mutant'

# Identify columns
Mutant_cols <- grep("^(TK57_|TK58_)", colnames(Mutant), value = TRUE)
Control_cols <- grep("^(TK43_|TK47_)", colnames(Mutant), value = TRUE)

# Sample TK43_ columns
if(length(Mutant_cols) >= 15282) {
  sampled_Mutant_cols <- sample(Mutant_cols, size = 15282, replace = FALSE)
} else {
  sampled_Mutant_cols<- Mutant_cols
  message("Not enough columns to sample, using all available.")
}

# Combine with all columns
final_columns <- c(sampled_Mutant_cols, Control_cols)

# Subset the data frame
Mutant_adjusted <- subset(Mutant, cells = final_columns)

save(Mutant_adjusted, file = "/media/thomaskim/Data/scATAC/Files/Dlx_Mutant_adjusted_Final.Robj")


####Plotting####

#Density
library(viridis)
library(viridisLite)
library(ggplot2)
Idents(Mutant_adjusted) <- "Cluster_Pass2"
my_levels <- c("Control","Mutant")
Grid <- Mutant_adjusted@reductions$umap@cell.embeddings
Grid <- data.frame(Grid)
Grid$orig.ident <- Mutant_adjusted@meta.data$Genotype
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

#Identify cluster
Idents(Mutant_adjusted) <- "Cluster_Pass2"
Cells1 <-WhichCells(Mutant_adjusted, idents = "AntID_ID_TT")
Cells1 <-WhichCells(Mutant_adjusted, idents = "SMN")
Cells1 <-WhichCells(Mutant_adjusted, idents = "MMN")
DimPlot(Mutant_adjusted, reduction = "umap", label = F, 
        pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
  ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf"))

#Percentage
library(dplyr)
library(reshape2)
library(Rmisc)
meta.data <- Mutant_adjusted@meta.data
Idents(Mutant_adjusted) <- "Cluster_Pass2"
graph <- 100*prop.table(table(Idents(Mutant_adjusted), Mutant_adjusted@meta.data$Genotype), margin = 1) #margin 1= row, 2 = column
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



TK43[['peakunion']] <- CreateAssayObject(counts = TK43.counts)
TK47[['peakunion']] <- CreateAssayObject(counts = TK47.counts)

DefaultAssay(TK43) <- "peakunion"
DefaultAssay(TK47) <- "peakunion"

