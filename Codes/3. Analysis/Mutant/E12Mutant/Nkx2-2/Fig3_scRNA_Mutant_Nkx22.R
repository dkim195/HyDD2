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
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram
setwd("/media/thomaskim/Data/")

####Nkx22_scRNA####
#plot
load(file = "scRNA/Robj/Nkx22_Final.Robj")
Idents(Mutant) <- "Cluster_Pass2"
Mutant <- RenameIdents(Mutant, "AntID_ID_TT" = "AntID_ID",
                       "LH (Pmch)" = "Tub",
                       "Prethalamus (Sp8)" = "PreThal",
                       "Prethalamus (Sp9)" = "PreThal",
                       "Prethalamus (Sst)" = "PreThal",
                       "PVH_SON (Cartpt)"= "PVH_SON", 
                       "PVH_SON (Otp)"= "PVH_SON",
                       "Tuberal (ARC_VMH)" = "Tub")
identities <- Idents(Mutant)
new_order <- sort(levels(identities))
Idents(Mutant) <- factor(identities, levels = new_order)
group <- c("AntID_ID" = "#808AA9",
           "Tub"= "#FF4136", "PreThal" = "#3A62D9", "PVH_SON"= "#DFCE93",
           "SMN"= "#93F6F8", "NPC" = "#F6B4F1",  
           "Neural Pro (Neurog2)" = "#4EF860",
           "Neural Pro" = "#F1CCB8", "Neural Pro (Ascl1)" = "#00FF1B",
           "PMN" = "#FFBB40", "MMN"= "#F0A881")

#UMAP
p <- DimPlot(Mutant, label = T, cols = group, raster = F)
png(filename = "Figures/Fig3/Nkx22/DEG_Nkx22_scRNA_raw.png", width = 900, height = 900)
print(p)
dev.off()

#UMAP
p <- DimPlot(Mutant, label = , cols = group, raster = F, split.by = "Genotype")+ NoLegend() + NoAxes()
png(filename = "Figures/Fig3/Nkx22/DEG_Nkx22_scRNA.png", width = 2400, height = 1200, res = 300)
print(p)
dev.off()

#split by mutant
check <- table(Mutant@active.ident)
check2 <- table(Mutant$Genotype)
Mutant <- AddMetaData(Mutant, Mutant@active.ident, "Cluster")
Mutant$CompositeLabel <- paste(Mutant$Cluster, Mutant$Genotype, sep = "_")
Idents(Mutant) <- "CompositeLabel"
new_order <-c("AntID_ID_Control", "AntID_ID_Mutant",
              "MMN_Control", "MMN_Mutant",
              "NPC_Control", "NPC_Mutant",
              "Neural Pro_Control", "Neural Pro_Mutant",
              "Neural Pro (Ascl1)_Control", "Neural Pro (Ascl1)_Mutant",
              "Neural Pro (Neurog2)_Control", "Neural Pro (Neurog2)_Mutant",
              "PMN_Control", "PMN_Mutant",
              "PVH_SON_Control", "PVH_SON_Mutant",
              "PreThal_Control", "PreThal_Mutant",
              "SMN_Control", "SMN_Mutant",
              "Tub_Control", "Tub_Mutant")
Idents(Mutant) <- factor(Idents(Mutant), levels = new_order)
markers <- FindAllMarkers(Mutant, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "Figures/Fig3/Nkx22/DEG_Nkx22_scRNA_Cluster_Genotype.csv")

#heatmap
library(pheatmap)
lists <- read.csv(file = "Figures/Fig3/Nkx22/DEG_Nkx22_scRNA_Cluster_Genotype.csv",
                  header = T, row.names = 1)
genes <- lists %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(avg_log2FC > 3) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  distinct(gene) %>%  # Ensure unique genes are selected if they appear in multiple clusters
  pull(gene)  # Extract the gene names as a character vector

#from DEG
genes <- c("Lhx8", "Prdm12", "D930028M14Rik", "Gabra2", "Pnoc", "Lhx6", 
           "B230323A14Rik", "Pcp4", "Foxb1", "Robo3", "Fgf13", "Vsnl1", 
           "Ednrb", "Pclaf", "Sgo1", "Gm32061", "Sfrp2", "Hes1", "Sgol1", 
           "Slc1a3", "Gsx2", "Helt", "Gsx1", "Gadd45g", "Plk3", "Neurog1", 
           "Rspo3", "Nhlh1", "Neurog2", "Sstr2", "Fn1", "Hmx3", "Hmx2", 
           "Cartpt", "Rxrg", "Crtac1", "Sst", "Sp8", "Meis2", "Sp9", "Gad2", 
           "C130021I20Rik", "Nr4a2", 
           "Irx3", "Irx5", "C1ql1", "Crnde", "Lmx1b", "Lmx1a", "Chchd10", 
           "Vgll2", "Nr5a1", "Nup62cl", "Sox14", "Nr0b1", "Pmch")
valid_genes <- genes[genes %in% rownames(Mutant[["SCT"]])]
Scaling  <- ScaleData(Mutant, features = valid_genes )
Avg.Scaling <- AverageExpression(Scaling, assay = "SCT", slot = "counts",
                                 verbose = T) #exponential minus 1 mean(expm1(x)) of the SCTTransform value
transposed_matrix <- t(Avg.Scaling$SCT[genes,])
scaled_transposed_matrix <- t(scale(t(transposed_matrix)))
newnames <- lapply(
  rownames(scaled_transposed_matrix),
  function(x) bquote(italic(.(x))))
heatmp <- pheatmap(scaled_transposed_matrix, fontsize = 8,scale = "column",
                   fontsize_col = 8, cluster_rows = FALSE,
                   cluster_cols = FALSE, cellwidth = 10, cellheight = 10, angle_col = 45,
                   color = colorRampPalette(c("blue", "white", "red"))(50),
                   labels_row = as.expression(newnames))
png(filename = "Figures/Fig3/Nkx22/heatmap_Nkx22_DEG_Italic.png", width = 3500, height = 1800, res = 300)
print(heatmp)
dev.off()
#newnames <- lapply(
#  rownames(Avg.Scaling$SCT[valid_genes,]),
#  function(x) bquote(italic(.(x))))
#
#p <- pheatmap(Avg.Scaling$SCT[valid_genes,], fontsize = 8, scale = "row",
#              fontsize_col = 8, cluster_rows = T,
#              cluster_cols = F, cellwidth = 10, cellheight = 10, angle_col = 45,
#              color = colorRampPalette(c("blue", "white", "red"))(50),
#              labels_row = as.expression(newnames))
#png(filename = "Figures/Fig3/Nkx22/heatmap_Nkx22_DEG_Italic.png", width = 3000, height = 2000)
#print(p)
#dev.off()
#p <- pheatmap(Avg.Scaling$SCT[valid_genes,], fontsize = 8, scale = "row",
#              fontsize_col = 8, cluster_rows = T,
#              cluster_cols = F, cellwidth = 10, cellheight = 10, angle_col = 45,
#              color = colorRampPalette(c("blue", "white", "red"))(50))
#png(filename = "Figures/Fig3/Nkx22/heatmap_Nkx22_DEG.png", width = 3000, height = 2000)
#print(p)
#dev.off()


#pattern DEG
genes <- c("Lhx8", "3110039M20Rik.1", "Lhx6", "Prdm12", "Col25a1", "Onecut3", 
           "Nkx6-1", "Nkx6-2", "Ppp1r14a", "Vsnl1", "Pax7", "Ntrk1", "Nppa", 
           "B230323A14Rik", "Foxb1", "Gadd45g", "Gsx2", "2010110K18Rik", 
           "Pif1", "Fam64a", "2810417H13Rik", "Hist2h3b", "Hist1h2bp", "Hist1h4n", 
           "Vtn", "Dkk2", "Lrtm1", "Msx1", "Pxdc1", "Halr1", "Wnt9a", "Calml4", 
           "Tbc1d4", "Aldh1l1", "Pmp22", "Has2", "Ptx3", "Arhgef6", "Col22a1", 
           "Tnc", "Hepacam", "Aldoc", "Tex15", "Gm29478", "Fgf8", "E030030I06Rik", 
           "Ctgf", "Tbx3", "Tfap2c", "Eomes", "Neurod4", "Mybpc1", "Neurog2", 
           "Sapcd2", "Gnat3", "Neurog1", "Htr3a", "BC049730", "Gm26771", 
           "Gm42047", "Dapk2", "Hist1h3c", "Hist1h2an", "Hist1h1b", "Hist1h2ap", 
           "Ufd1l", "Erdr1", "Mcam", "Atp1a2", "Hey2", "Gsx1", "Rgs16", 
           "Helt", "Gm11549", "Msx3", "Nbl1", "Fli1", "Optc", "Sstr2", "Actc1", 
           "Palmd", "Ebf1", "Hmx2", "Hmx3", "Rrad", "Gal", "Gm48336", "Naip2", 
           "Ccrl2", "Gstp2", "Dpcr1", "Fam101a", "Cartpt", "Gbx2", "Sema3c", 
           "Nfam1", "Car8", "Crabp1", "Sncg", "Pou4f1", "Otp", "Rxrg", "Chac1", 
           "Insm2", "Calcr", "Unc13c", "Sphkap", "Slc30a3", "Sp9", "Sp8", 
           "2700033N17Rik", "Meis2", "G630016G05Rik", "Foxa2", "Fignl2", 
           "AC119932.1", "Irx6", "Bfsp2", "C130021I20Rik", "Gm26872", "Lmx1b", 
           "Calb2", "Nr4a2", "Neurod6", "Sst", "Ret", "Pcdh20", "Npy", "Shox2", 
           "Lhx9", "Prdm13", "Sox14", "Nup62cl", "Vgll2", "Nr5a1", "Thsd4", 
           "Gpr149", "Cfap77", "Wnt8b", "Sp5", "Bmp4", "Tagln2", "Maob", 
           "Mboat1")
valid_genes <- genes[genes %in% rownames(Mutant[["SCT"]])]
Scaling  <- ScaleData(Mutant, features = valid_genes )
Avg.Scaling <- AverageExpression(Scaling, assay = "SCT", slot = "counts",
                                 verbose = T) #exponential minus 1 mean(expm1(x)) of the SCTTransform value
transposed_matrix <- t(Avg.Scaling$SCT[valid_genes,])
scaled_transposed_matrix <- t(scale(t(transposed_matrix)))
newnames <- lapply(
  rownames(scaled_transposed_matrix),
  function(x) bquote(italic(.(x))))
heatmp <- pheatmap(scaled_transposed_matrix, fontsize = 8,scale = "column",
                   fontsize_col = 8, cluster_rows = FALSE,
                   cluster_cols = FALSE, cellwidth = 10, cellheight = 10, angle_col = 45,
                   color = colorRampPalette(c("blue", "white", "red"))(50),
                   labels_row = as.expression(newnames))
png(filename = "Figures/Fig3/Nkx22/heatmap_Nkx22_PatternGene_Italic.png", width = 8000, height = 1800, res = 300)
print(heatmp)
dev.off()

#newnames <- lapply(
#  rownames(Avg.Scaling$SCT[valid_genes,]),
#  function(x) bquote(italic(.(x))))

#p <- pheatmap(Avg.Scaling$SCT[valid_genes,], fontsize = 8, scale = "row",
#              fontsize_col = 8, cluster_rows = T,
#              cluster_cols = F, cellwidth = 10, cellheight = 10, angle_col = 45,
#              color = colorRampPalette(c("blue", "white", "red"))(50),
#              labels_row = as.expression(newnames))
#png(filename = "Figures/Fig3/Nkx22/heatmap_Nkx22_PatternGene_Italic.png", width = 3000, height = 2000)
#print(p)
#dev.off()

#p <- pheatmap(Avg.Scaling$SCT[valid_genes,], fontsize = 8, scale = "row",
#              fontsize_col = 8, cluster_rows = T,
#              cluster_cols = F, cellwidth = 10, cellheight = 10, angle_col = 45,
#              color = colorRampPalette(c("blue", "white", "red"))(50))
#png(filename = "Figures/Fig3/Nkx22/heatmap_Nkx22_PatternGene.png", width = 3000, height = 2000)
#print(p)
#dev.off()

#distribution
library(viridis)
library(viridisLite)
library(ggplot2)
library(cetcolor)
Idents(Mutant) <- "Cluster_Pass2"
Mutant <- RenameIdents(Mutant, "AntID_ID_TT" = "AntID_ID",
                       "LH (Pmch)" = "Tub",
                       "Prethalamus (Sp8)" = "PreThal",
                       "Prethalamus (Sp9)" = "PreThal",
                       "Prethalamus (Sst)" = "PreThal",
                       "PVH_SON (Cartpt)"= "PVH_SON", 
                       "PVH_SON (Otp)"= "PVH_SON",
                       "Tuberal (ARC_VMH)" = "Tub")
# Generate UMAP plot
pl1 <- DimPlot(Mutant, combine = F, cols = group, split.by = "Genotype")
# Custom color scale
scale.col <- cet_pal(16, name = "fire")
# Function to threshold density values
threshold_density <- function(x, from = 0, to = 100) {
  x <- scales::rescale(x, to = c(from, to))
  x[x < from] <- from
  x[x > to] <- to
  return(x)
}
# Make plot with thresholded density values
gplot <- pl1[[1]] & 
  stat_density_2d(aes_string(x = "UMAP_1", y = "UMAP_2", fill = "threshold_density(after_stat(level))"), 
                  linewidth = 0.2, geom = "density_2d_filled", 
                  colour = "ivory", alpha = 0.4, n = 150, h = c(1.2, 1.2)) & 
  scale_fill_gradientn(colours = scale.col, limits = c(0, 100)) & 
  DarkTheme()
ggsave(filename = "Figures/Fig3/Nkx22/densityheatmap_Nkx22.eps",
       plot = gplot, width = 20, height = 10, units = "in", device = cairo_ps)

#identities <- Idents(Mutant)
#Grid <- Mutant@reductions$umap@cell.embeddings
#Grid <- data.frame(Grid)
#Grid$orig.ident <- Mutant$Genotype
#Grid_ds <- Grid[,1:2]
#p <- ggplot(Grid, aes(x = UMAP_1, y = UMAP_2)) +
#  geom_point(data=Grid_ds, size=0.1, alpha=0.1, color="white") +
#  scale_fill_viridis(option="A", name = "Density")+
#  facet_grid(~orig.ident) +
#  stat_density_2d(geom="raster",aes(fill=stat(ndensity)), contour=F) +
#  theme_classic() +
#  scale_x_continuous(expand=c(0,0)) +
#  scale_y_continuous(expand=c(0,0)) +
#  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 16),
#        strip.text.x = element_text(size = 24)) +
#  coord_fixed()
#png(filename = "Figures/Fig3/Nkx22/densityheatmap_Nkx22.png", width = 1200, height = 1000)
#print(p)
#dev.off()

#density
library(dplyr)
library(reshape2)
library(Rmisc)
Idents(Mutant) <- "Cluster_Pass2"
Mutant <- RenameIdents(Mutant, "AntID_ID_TT" = "AntID_ID",
                       "LH (Pmch)" = "Tub",
                       "Prethalamus (Sp8)" = "PreThal",
                       "Prethalamus (Sp9)" = "PreThal",
                       "Prethalamus (Sst)" = "PreThal",
                       "PVH_SON (Cartpt)"= "PVH_SON", 
                       "PVH_SON (Otp)"= "PVH_SON",
                       "Tuberal (ARC_VMH)" = "Tub")
identities <- Idents(Mutant)
new_order <- sort(levels(identities))
Idents(Mutant) <- factor(identities, levels = new_order)
meta.data <- Mutant@meta.data
graph <- 100*prop.table(table(Idents(Mutant), Mutant$Genotype), margin = 1) #margin 1= row, 2 = column
dat <- melt(graph)
colors <- brewer.pal(8, "Set2")
# Your ggplot code with the new colors
gplot <- ggplot(data = dat, aes(x = factor(Var1), y = value, fill = Var2)) + 
  geom_bar(stat = "identity", colour = "black", width = 1) +
  scale_fill_manual(values = colors) +
  ylab("Cluster proportion") + xlab("Cluster identity") +
  ggtitle("Proportion of cell types") +
  guides(fill = guide_legend(title = "Group")) +
  theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())
ggsave(filename = "Figures/Fig3/Nkx22/densitypercentage_Nkx22.eps",
       plot = gplot, width = 20, height = 10, units = "in", device = cairo_ps)
#png(filename = "Figures/Fig3/Nkx22/densitypercentage_Nkx22.png", width = 1200, height = 1000)
#print(p)
#dev.off()



####identify cluster####
Cells1 <-WhichCells(Mutant, idents = "AntID_ID_TT")
DimPlot(Mutant, reduction = "umap", label = F, 
        pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
  ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf"))


