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
load(file = "scATAC/Files/TK53.Robj")
Idents(TK53) <- "orig.ident"
TK53 <- RenameIdents(TK53, "ATAC" = "E16")
TK53 <- AddMetaData(TK53, TK53@active.ident, "Age")
Idents(TK53) <- "orig.ident"
TK53 <- RenameIdents(TK53, "ATAC" = "TK53")
TK53 <- AddMetaData(TK53, TK53@active.ident, "Sample")
Idents(TK53) <- "orig.ident"
TK53 <- RenameIdents(TK53, "ATAC" = "Control")
TK53 <- AddMetaData(TK53, TK53@active.ident, "Genotype")
Idents(TK53) <- "Age"

load(file = "scATAC/Files/TK54.Robj")
Idents(TK54) <- "orig.ident"
TK54 <- RenameIdents(TK54, "ATAC" = "E16")
TK54 <- AddMetaData(TK54, TK54@active.ident, "Age")
Idents(TK54) <- "orig.ident"
TK54 <- RenameIdents(TK54, "ATAC" = "TK54")
TK54 <- AddMetaData(TK54, TK54@active.ident, "Sample")
Idents(TK54) <- "orig.ident"
TK54 <- RenameIdents(TK54, "ATAC" = "Control")
TK54 <- AddMetaData(TK54, TK54@active.ident, "Genotype")
Idents(TK54) <- "Age"

#mutant
load(file = "scATAC/Files/TK55.Robj")
Idents(TK55) <- "orig.ident"
TK55 <- RenameIdents(TK55, "ATAC" = "E16")
TK55 <- AddMetaData(TK55, TK55@active.ident, "Age")
Idents(TK55) <- "orig.ident"
TK55 <- RenameIdents(TK55, "ATAC" = "TK55")
TK55 <- AddMetaData(TK55, TK55@active.ident, "Sample")
Idents(TK55) <- "orig.ident"
TK55 <- RenameIdents(TK55, "ATAC" = "Mutant")
TK55 <- AddMetaData(TK55, TK55@active.ident, "Genotype")
Idents(TK55) <- "Age"

load(file = "scATAC/Files/TK56.Robj")
Idents(TK56) <- "orig.ident"
TK56 <- RenameIdents(TK56, "ATAC" = "E16")
TK56 <- AddMetaData(TK56, TK56@active.ident, "Age")
Idents(TK56) <- "orig.ident"
TK56 <- RenameIdents(TK56, "ATAC" = "TK56")
TK56 <- AddMetaData(TK56, TK56@active.ident, "Sample")
Idents(TK56) <- "orig.ident"
TK56 <- RenameIdents(TK56, "ATAC" = "Mutant")
TK56 <- AddMetaData(TK56, TK56@active.ident, "Genotype")
Idents(TK56) <- "Age"

####Fragment####
DefaultAssay(TK53) <- 'peaks'
DefaultAssay(TK54) <- 'peaks'
DefaultAssay(TK55) <- 'peaks'
DefaultAssay(TK56) <- 'peaks'

frags <- Fragments(TK53)  # get list of fragment objects
Fragments(TK53) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK53/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK53) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK54)  # get list of fragment objects
Fragments(TK54) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK54/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK54) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK55)  # get list of fragment objects
Fragments(TK55) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK55/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK55) <- frags  # assign update list of fragment objects back to the assay

frags <- Fragments(TK56)  # get list of fragment objects
Fragments(TK56) <- NULL  # remove fragment information from assay
newpath <- "scATAC/Fragments/TK56/fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(TK56) <- frags  # assign update list of fragment objects back to the assay

####Intersection####
# Find the union of peaks between the two datasets
# Extract peak regions from both datasets
DefaultAssay(TK53) <- 'peaks'
DefaultAssay(TK54) <- 'peaks'
DefaultAssay(TK55) <- 'peaks'
DefaultAssay(TK56) <- 'peaks'

peaks4 <- GetAssayData(TK53, assay = "peaks", slot = "ranges")
peaks2 <- GetAssayData(TK54, assay = "peaks", slot = "ranges")
peaks3 <- GetAssayData(TK55, assay = "peaks", slot = "ranges")
peaks1 <- GetAssayData(TK56, assay = "peaks", slot = "ranges")

# Find the intersection of peaks
combined.peaks <- reduce(x = c(peaks1, peaks2, peaks3, peaks4))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# Subset both objects to include only common peaks
DefaultAssay(TK53) <- 'peaks'
DefaultAssay(TK54) <- 'peaks'
DefaultAssay(TK55) <- 'peaks'
DefaultAssay(TK56) <- 'peaks'

TK53.counts <- FeatureMatrix(
  fragments = Fragments(TK53),
  features = combined.peaks,
  cells = colnames(TK53))
TK54.counts <- FeatureMatrix(
  fragments = Fragments(TK54),
  features = combined.peaks,
  cells = colnames(TK54))
TK55.counts <- FeatureMatrix(
  fragments = Fragments(TK55),
  features = combined.peaks,
  cells = colnames(TK55))
TK56.counts <- FeatureMatrix(
  fragments = Fragments(TK56),
  features = combined.peaks,
  cells = colnames(TK56))

TK53[['peakunion']] <- CreateAssayObject(counts = TK53.counts)
TK54[['peakunion']] <- CreateAssayObject(counts = TK54.counts)
TK55[['peakunion']] <- CreateAssayObject(counts = TK55.counts)
TK56[['peakunion']] <- CreateAssayObject(counts = TK56.counts)

DefaultAssay(TK56) <- 'peakunion'
DefaultAssay(TK54) <- 'peakunion'
DefaultAssay(TK55) <- 'peakunion'
DefaultAssay(TK53) <- 'peakunion'

rm(TK56.counts)
rm(TK54.counts)
rm(TK55.counts)
rm(TK53.counts)

####Merging####
# Merge the datasets
Mutant <- merge(TK56, y = c(TK54, TK55, TK53),
                add.cell.ids = c("TK56", "TK54", "TK55", "TK53"))

rm(TK56)
rm(TK54)
rm(TK55)
rm(TK53)

save(Mutant, file = "/media/thomaskim/Data/scATAC/Files/Foxd1_Mutant_adjusted.Robj")

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
DimPlot(Mutant, label = , group.by = "Genotype")

####Transfer Label####
#RNA
load(file = "/media/thomaskim/Data/scRNA/Robj/E16_Atlas.Robj")
E16 <- UpdateSeuratObject(E16)
E16@meta.data$old.ident <- NULL
E16 <- AddMetaData(E16, E16@active.ident, "Cluster_Pass2")
E16@meta.data$Cluster <- NULL
E16@meta.data$SCT_snn_res.2 <- NULL
E16@meta.data$seurat_clusters <- NULL
Idents(E16) <- "Age"

E16 <- RenameIdents(E16, "E16_Rep1" = "Atlas", "E16_Rep2" = "Atlas")
E16 <- AddMetaData(E16, E16@active.ident, "Reference")
Idents(E16) <- "Cluster_Pass2"

E16 <- RenameIdents(E16, "Avp_Gal_Oxt (PVN_SON)" = "PVN_SON", "Bsx (PMN)"= "PMN", 
                    "Calb2 (SMN)" = "SMN", "Cck (MMN)" = "MMN",
                    "Check (POA_SCN)" = "POA_SCN", "Ghrh_Gal (PMN)"= "PMN", 
                    "Npvf_Grp_Hcrt (LH)" = "LH", "Npy_Pomc_Agrp (ARC)" = "ARC_VMH",
                    "NPC (Astro)" = "NPC (Glial)", 
                    "NPC (Epen)" = "NPC (Glial)", "NPC (Glial)" = "NPC (Glial)",
                    "NPC (Oligo)" = "NPC (Glial)", "NPC (Tanycyte)" = "NPC (Glial)",
                    "Nts_Tac1_Cck (MMN)" = "MMN", 
                    "Penk (POA_SCN)"= "POA_SCN", "Pmch_Trh (LH)"= "LH",
                    "Rora_Rorb (POA_SCN)"= "POA_SCN", "Tac2 (ARC_VMH)"= "ARC_VMH")

E16 <- AddMetaData(E16, E16@active.ident, "Cluster_Pass2")
Idents(E16) <- "Cluster_Pass2"

DefaultAssay(Mutant) <- 'RNA'
DefaultAssay(E16) <- 'RNA'

E16 <- FindVariableFeatures(object = E16,  nfeatures = 350)

transfer.anchors <- FindTransferAnchors(
  reference = E16, query = Mutant,
  normalization.method = "LogNormalize",
  reduction = "pcaproject", dims = 1:30)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = E16$Cluster_Pass2,
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

rm(E16)
save(Mutant, file = "/media/thomaskim/Data/scATAC/Files/Foxd1_Mutant_adjusted.Robj")

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

####Match cell number####
Mutant$Cluster_Pass2 <- Mutant$predicted.id
Idents(Mutant) <- "Genotype"

# Assuming 'Mutant' is your Seurat object
# Assume you have a metadata column identifying each cell as 'control' or 'mutant'

# Identify columns
Mutant_cols <- grep("^(TK55_|TK56_)", colnames(Mutant), value = TRUE)
Control_cols <- grep("^(TK53_|TK54_)", colnames(Mutant), value = TRUE)

# Sample TK56_ columns
if(length(Mutant_cols) >= 15000) {
  sampled_Mutant_cols <- sample(Mutant_cols, size = 15000, replace = FALSE)
} else {
  sampled_Mutant_cols<- Mutant_cols
  message("Not enough columns to sample, using all available.")
}

# Combine with all columns
final_columns <- c(sampled_Mutant_cols, Control_cols)

# Subset the data frame
Mutant_adjusted <- subset(Mutant, cells = final_columns)
table(Mutant_adjusted$Genotype)
rm(Mutant)
save(Mutant_adjusted, file = "/media/thomaskim/Data/scATAC/Files/Foxd1_Mutant_adjusted_Final.Robj")


library(dplyr)
library(reshape2)
library(Rmisc)
meta.data <- Mutant_adjusted@meta.data
Idents(Mutant_adjusted) <- "predicted.id"
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



TK56[['peakunion']] <- CreateAssayObject(counts = TK56.counts)
TK54[['peakunion']] <- CreateAssayObject(counts = TK54.counts)

DefaultAssay(TK56) <- "peakunion"
DefaultAssay(TK54) <- "peakunion"

