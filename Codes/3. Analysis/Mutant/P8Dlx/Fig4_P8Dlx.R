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
library(cetcolor)
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram
setwd("/media/thomaskim/Data/")

####Dlx_scRNA####
#plot
load(file = "scRNA/Robj/P8_Dlx_Final.Robj")

markers <- FindAllMarkers(Mutant, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "Figures/Fig4/P8_Dlx_DEG.csv")

DimPlot(Mutant, label = T, group.by = "Cluster_Pass3") + NoLegend()

Idents(Mutant) <- "Cluster_Pass3"
#modifying based on 1. location and neuropeptide
#Annotation from Mutant P8 is not super great
Mutant <- RenameIdents(Mutant, "ARC (Npw_Agrp)" = "ARC (Npy_Agrp)", 
                       "ARC (Npy_Agrp)" = "ARC (Npy_Agrp)", 
                       "ARC (Pomc)" = "ARC (Pomc)", 
                       "DMH (Npw_PMN)" = "DMH_PMN (Npw)",
                       "DMH-LH (Grp_Cck)" = "DMH_PMN (Grp_Bsx)", 
                       "LH (Hcrt)" = "LH_DMH (Hcrt_Npvf)",
                       "LH (Pmch_Tac2)" = "ARC_LH (Kiss1_Tac2_Pmch)", 
                       "LH (Pmch_Trh)" = "PVH_SON (Oxt_Trh_Crh)",
                       "MMN (Calb1_Cck)" = "MMN", 
                       "MMN (Cartpt_Cck)" = "MMN" ,
                       "MMN (Cck)" = "MMN", 
                       "MMN (Nts_Cartpt)" = "MMN",
                       "MMN (Nts_Tac1_Cck)" = "MMN", 
                       "MMN (Pcp4_Cck)" = "MMN",
                       "PMN (Ghrh_Gal_Oxt)" = "ARC_PMN (Ghrh_Gal_Th)" ,
                       "PreThal" = "PreThal",
                       "PVN_SON (Avp_Gal_Oxt)" = "PVN_SON (Avp)", 
                       "SCN (Avp)" = "PVN_SON (Avp)", 
                       "SCN (Rorb)" = "SCN (Nms)",
                       "SCN (Vip)" = "SCN (Vip)",
                       "SMN (Calb2_Tac1)" = "SMN",
                       "SMN (Calb2)" = "Thalamus", 
                       "SMN (Nts)" = "SMN",
                       "Sst (ARC-TMN?)" = "UNKNOWN (Sst_Pthlh)" ,
                       "Sst (TMN-PH?)" = "TMN_PH (Sst_Thrb)",
                       "Sst_Npy"  = "ARC (Sst_Npy_Th)",
                       "Tac2" = "SMN",
                       "TMN_PH (Hdc)" = "TMN_PH (Hdc_Prph)",
                       "VMH"  = "VMH (Nr5a1)",
                       "VMH (Tac1)"  = "VMH (Nr5a1_Tac1)")
identities <- Idents(Mutant)
new_order <- sort(levels(identities))
Idents(Mutant) <- factor(identities, levels = new_order)
group <- c(
  "ARC (Npy_Agrp)" = "#FB8072", 
  "ARC (Pomc)" = "#B3DE69",
  "ARC (Sst_Npy_Th)" = "#CE6ACB", 
  "ARC_LH (Kiss1_Tac2_Pmch)" = "#F09C6B",
  "ARC_PMN (Ghrh_Gal_Th)" = "#E78AC3",
  "DMH_PMN (Grp_Bsx)"= "#BEBADA",
  "DMH_PMN (Npw)" = "#FFFFB3", 
  "LH_DMH (Hcrt_Npvf)"= "#DD8F00",
  "MMN"= "#984EA3",
  "PreThal" = "#E5C494",
  "PVH_SON (Oxt_Trh_Crh)" = "#1F78B4",
  "PVN_SON (Avp)" = "#E9DEC1",
  "SCN (Nms)" = "#33A02C",
  "SCN (Vip)" = "#AA102E",
  "SMN"= "#FB9A99",
  "Thalamus" = "#CAB2D6",
  "TMN_PH (Hdc_Prph)"= "#6A3D9A",
  "TMN_PH (Sst_Thrb)"= "#FF7F00",
  "UNKNOWN (Sst_Pthlh)"= "#BC80BD", 
  "VMH (Nr5a1_Tac1)"= "#FCCDE5",
  "VMH (Nr5a1)" = "#FDB462"
)

#UMAP
p <- DimPlot(Mutant, label = T, cols = group, raster = F)
png(filename = "Figures/Fig4/P8_Dlx_scRNA_raw.png", width = 900, height = 900)
print(p)
dev.off()

#UMAP
p <- DimPlot(Mutant, label = , cols = group, raster = F, split.by = "Genotype")+ NoLegend() + NoAxes()
png(filename = "Figures/Fig4/P8_Dlx_scRNA.png", width = 1200, height = 1000)
print(p)
dev.off()

p <- DimPlot(Mutant, label = , cols = group, raster = F, split.by = "Genotype")+ NoLegend() + NoAxes()
png(filename = "Figures/Fig4/P8_Dlx_scRNA2.png", width = 2400, height = 1200, res = 300)
print(p)
dev.off()


#split by mutant
check <- table(Mutant@active.ident)
check2 <- table(Mutant$Genotype)
Mutant <- AddMetaData(Mutant, Mutant@active.ident, "Cluster")
Mutant$CompositeLabel <- paste(Mutant$Cluster, Mutant$Genotype, sep = "_")
Idents(Mutant) <- "CompositeLabel"

clusters <- as.data.frame(table(Mutant$Cluster))
genotypes <- as.data.frame(table(Mutant$Genotype))
all_combinations <- expand.grid(Cluster = clusters$Var1, Genotype = genotypes$Var1)
all_combinations$Combined <- paste(all_combinations$Cluster, all_combinations$Genotype, sep = "_")
all_combinations <- all_combinations[order(all_combinations$Cluster, all_combinations$Genotype), ]
new_order <- all_combinations$Combined
new_order
Idents(Mutant) <- factor(Idents(Mutant), levels = new_order)
markers <- FindAllMarkers(Mutant, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "Figures/Fig4/P8_Dlx_scRNA_Cluster_Genotype.csv")

#heatmap
library(pheatmap)
lists <- read.csv(file = "Figures/Fig4/P8_Dlx_scRNA_Cluster_Genotype.csv",
                  header = T, row.names = 1)
genes <- lists %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(avg_log2FC > 3) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  distinct(gene) %>%  # Ensure unique genes are selected if they appear in multiple clusters
  pull(gene)  # Extract the gene names as a character vector

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
png(filename = "Figures/Fig4/P8_heatmap_Dlx_DEG_Italic.png", width = 5000, height = 2500, res = 300)
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
#png(filename = "Figures/Fig4/P8_heatmap_Dlx_DEG_Italic.png", width = 3000, height = 2000)
#print(p)
#dev.off()

#p <- pheatmap(Avg.Scaling$SCT[valid_genes,], fontsize = 8, scale = "row",
#              fontsize_col = 8, cluster_rows = T,
#              cluster_cols = F, cellwidth = 10, cellheight = 10, angle_col = 45,
#              color = colorRampPalette(c("blue", "white", "red"))(50))
#png(filename = "Figures/Fig4/P8_heatmap_Dlx_DEG.png", width = 3000, height = 2000)
#print(p)
#dev.off()

#pattern DEG
library(pheatmap)
lists <- read.csv(file = "Figures/Fig1/Pattern_scRNA//DEG_pattern_scRNA.csv",
                  header = T, row.names = 1)
#selected groups
select_group <- c("AntID_ID", "MMN", "MMN (Prog)", "Neural Pro (Ascl1)", 
                  "Neural Pro (Neurog2)","PMN", 
                  "PMN (Gsx1)", "PMN (Pro)", "POA/SCN", "PreThal (Pro)", "PreThal (Sp8)", 
                  "PreThal (Sp9)", "PreThal_AntID", "PVH_SON", "PVH_SON_AntVen", 
                  "SMN", "SMN (Prog)", "Tub (ARC_VMH)", "Tub (LH-Pmch)", 
                  "Tub (LH)", "Tub (Prog)")
genes <- lists %>%
  filter(cluster %in% select_group) %>%  # Keep only select_group
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(avg_log2FC > 3) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  distinct(gene) %>%  # Ensure unique genes are selected if they appear in multiple clusters
  pull(gene)  # Extract the gene names as a character vector
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
png(filename = "Figures/Fig4/P8_heatmap_Dlx_PatternGene_Italic.png", width = 3500, height = 2500, res = 300)
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
#png(filename = "Figures/Fig4/P8_heatmap_Dlx_PatternGene_Italic.png", width = 3000, height = 2000)
#print(p)
#dev.off()

#p <- pheatmap(Avg.Scaling$SCT[valid_genes,], fontsize = 8, scale = "row",
#              fontsize_col = 8, cluster_rows = T,
#              cluster_cols = F, cellwidth = 10, cellheight = 10, angle_col = 45,
#              color = colorRampPalette(c("blue", "white", "red"))(50))
#png(filename = "Figures/Fig4/P8_heatmap_Dlx_PatternGene.png", width = 3000, height = 2000)
#print(p)
#dev.off()

#neurogenesis DEG
library(pheatmap)
lists <- read.csv(file = "Figures//Fig1/Neurogenesis_scRNA//DEG_neurogenesis_scRNA_compact.csv",
                  header = T, row.names = 1)
#all groups
genes <- lists %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(avg_log2FC > 3) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  distinct(gene) %>%  # Ensure unique genes are selected if they appear in multiple clusters
  pull(gene)  # Extract the gene names as a character vector
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
png(filename = "Figures/Fig4/P8_heatmap_Dlx_NeurogenesisGene_Italic.png", width = 3500, height = 2500, res = 300)
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
#png(filename = "Figures/Fig4/P8_heatmap_Dlx_NeurogenesisGene_Italic.png", width = 3000, height = 2000)
#print(p)
#dev.off()

#p <- pheatmap(Avg.Scaling$SCT[valid_genes,], fontsize = 8, scale = "row",
#              fontsize_col = 8, cluster_rows = T,
#              cluster_cols = F, cellwidth = 10, cellheight = 10, angle_col = 45,
#              color = colorRampPalette(c("blue", "white", "red"))(50))
#png(filename = "Figures/Fig4/P8_heatmap_Dlx_NeurogenesisGene.png", width = 3000, height = 2000)
#print(p)
#dev.off()

#distribution
library(viridis)
library(viridisLite)
library(ggplot2)
#Idents(Mutant) <- "Cluster"
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
#png(filename = "Figures/Fig4/P8_densityheatmap_Dlx.png", width = 1200, height = 1000)
#print(p)
#dev.off()

#distribution 2
library(cetcolor)
library(Seurat)
library(ggplot2)
# generate UMAP plot
#thresholding
# Generate UMAP plot
pl1 <- DimPlot(Mutant, combine = F, cols = group, split.by = "Genotype")
pl1 <- lapply(pl1, function(p) p + theme(legend.position = "none"))
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
ggsave(filename = "Figures/Fig4/P8_densityheatmap_Dlx2.eps",
       plot = gplot, width = 20, height = 10, units = "in", device = cairo_ps)

#pl2 <- pl1[[1]] & 
#  stat_density_2d(aes_string(x = "UMAP_1", y = "UMAP_2", fill = "threshold_density(after_stat(level))"), 
#                  linewidth = 0.2, geom = "density_2d_filled", 
#                  colour = "ivory", alpha = 0.4, n = 150, h = c(1.2, 1.2)) & 
#  scale_fill_gradientn(colours = scale.col, limits = c(0, 100)) & 
#  DarkTheme()
#png(filename = "Figures/Fig4/P8_densityheatmap_Dlx2.png", width = 1200, height = 1000)
#print(pl2)
#dev.off()

#density
library(dplyr)
library(reshape2)
library(Rmisc)
#Idents(Mutant) <- "Cluster"
#identities <- Idents(Mutant)
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
ggsave(filename = "Figures/Fig4/P8_densitypercentage_Dlx.eps",
       plot = gplot, width = 20, height = 10, units = "in", device = cairo_ps)

#png(filename = "Figures/Fig4/P8_densitypercentage_Dlx.png", width = 1200, height = 1000)
#print(p)
#dev.off()



####identify cluster####
Cells1 <-WhichCells(Mutant, idents = "AntID_ID_TT")
DimPlot(Mutant, reduction = "umap", label = F, 
        pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
  ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf"))


