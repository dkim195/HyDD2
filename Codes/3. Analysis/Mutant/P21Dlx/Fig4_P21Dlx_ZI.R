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

####P21_Dlx_sNRNA####
#plot
load(file = "scRNA/Robj/DLXCKO_P21_ZITRN_Neuron_Final.Robj")

markers <- FindAllMarkers(DLXCKO, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "Figures/Fig5/P21_ZI_Dlx_DEG.csv")

DimPlot(DLXCKO, label = T, group.by = "Clusters") + NoLegend()

Idents(DLXCKO) <- "Clusters"
#modifying based on 1. location and neuropeptide
#Annotation from DLXCKO P21_ZI is not super great
identities <- Idents(DLXCKO)
new_order <- sort(levels(identities))
Idents(DLXCKO) <- factor(identities, levels = new_order)
group <- c("Glutamatergic"= "#33A02C",
           "Lef1+"= "#FF7F00",
           "Calb1+"= "#B3DE69",
           "Sst+"= "#FFFFB3", 
           "Meis2+" = "#1F78B4",
           "Lhx1os+"= "#FCCDE5", 
           "Nfib+"= "#FF4136", 
           "Dlx6os1+" = "#3A62D9",
           "Pax6+"= "#93F6F8",
           "Oxt+"= "#984EA3",
           "Hcrt+"= "#DD8F00",
           "Crh+" = "#E9DEC1",
           "Foxg1+" = "#CAB2D6",
           "Tbr1+"= "#BC80BD", 
           "Pmch+"  = "#00FF1B")

#UMAP
p <- DimPlot(DLXCKO, label = T, cols = group, raster = F)
png(filename = "Figures/Fig5/P21_ZI_Dlx_scRNA_raw.png", width = 900, height = 900)
print(p)
dev.off()

#UMAP
p <- DimPlot(DLXCKO, label = , cols = group, raster = F, split.by = "Genotype_Final", pt.size = 2)+ NoLegend() + NoAxes()
png(filename = "Figures/Fig5/P21_ZI_Dlx_scRNA.png", width = 1200, height = 1000)
print(p)
dev.off()

p <- DimPlot(DLXCKO, label = , cols = group, raster = F, split.by = "Genotype_Final", pt.size = 0.5)+ NoLegend() + NoAxes()
png(filename = "Figures/Fig5/P21_ZI_Dlx_scRNA2.png", width = 2400, height = 1200, res = 300)
print(p)
dev.off()


#split by DLXCKO
check <- table(DLXCKO@active.ident)
check2 <- table(DLXCKO$Genotype_Final)
DLXCKO <- AddMetaData(DLXCKO, DLXCKO@active.ident, "Cluster")
DLXCKO$CompositeLabel <- paste(DLXCKO$Cluster, DLXCKO$Genotype_Final, sep = "_")
Idents(DLXCKO) <- "CompositeLabel"

clusters <- as.data.frame(table(DLXCKO$Cluster))
genotypes <- as.data.frame(table(DLXCKO$Genotype_Final))
all_combinations <- expand.grid(Cluster = clusters$Var1, Genotype = genotypes$Var1)
all_combinations$Combined <- paste(all_combinations$Cluster, all_combinations$Genotype, sep = "_")
all_combinations <- all_combinations[order(all_combinations$Cluster, all_combinations$Genotype), ]
new_order <- all_combinations$Combined
new_order
Idents(DLXCKO) <- factor(Idents(DLXCKO), levels = new_order)
markers <- FindAllMarkers(DLXCKO, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "Figures/Fig5/P21_ZI_Dlx_scRNA_Cluster_Genotype.csv")

#heatmap
library(pheatmap)
lists <- read.csv(file = "Figures/Fig5/P21_ZI_Dlx_scRNA_Cluster_Genotype.csv",
                  header = T, row.names = 1)
genes <- lists %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(avg_log2FC > 3) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  distinct(gene) %>%  # Ensure unique genes are selected if they appear in multiple clusters
  pull(gene)  # Extract the gene names as a character vector

genes <- c("C6", "Gxylt2", "Ntf3", "Gm29478", "Cfh", "Ankrd31", "Adgrg6", 
           "Bche", "Maob", "Pla2r1", "Avp", "Agtr1a", "Fosb", "Uncx", "Acot12", 
           "Cxcl2", "Cxcl1", "Crh", "Gm38505", "Sox2ot", "Gm42695", "Gm45680", 
           "Grp", "Gpr179", "Gm50024", "Sh3rf2", "Adora2a", "6530403H02Rik", 
           "Rarb", "Crym", "9130024F11Rik", "Foxg1", "1700003D09Rik", "Npvf", 
           "Slc25a48", "4930548K13Rik", "Csta2", "Rfx2", "Gm30302", "Vgll2", 
           "Hcrt", "Krt20", "Gm35696", "Car4", "Gm14866", "Adamts19", "Wfdc10", 
           "Lrrc63", "Pdlim1", "Greb1", "Crispld2", "Gna14", "Gm20713", 
           "Gpc5", "4930505G20Rik", "Gata3", "Gm15339", "Colq", "Eepd1", 
           "Lrrc38", "Gm2115", "Crtam", "Sapcd1", "Gm12153", "Adamtsl5", 
           "Gm14862", "Otx2os1", "E130114P18Rik", "Otx2", "Gm534", "Gm20098", 
           "Gm50423", "Oxt", "Esr2", "Gm10475", "Sim1", "Gm33203", "Gm40663", 
           "Gm36723", "Gm48821", "Pdgfd", "Pax6", "4930554C24Rik", "Gm14015", 
           "Gm28342", "Otx1", "Pmch", "Tacr3", "Cartpt", "Mup6", "Gm31804", 
           "A2m", "Corin", "Gm20647", "Slc9a2", "Samd3", "Cngb3", "8430419K02Rik", 
           "Gm12023", "Gm16141","Gal","Atoh7","Pnoc","Penk","Pvalb")


valid_genes <- genes[genes %in% rownames(DLXCKO[["SCT"]])]
Scaling  <- ScaleData(DLXCKO, features = valid_genes )
Avg.Scaling <- AverageExpression(Scaling, assay = "SCT", slot = "counts",
                                 verbose = T) #exponential minus 1 mean(expm1(x)) of the SCTTransform value

transposed_matrix <- t(Avg.Scaling$SCT[valid_genes,])
scaled_transposed_matrix <- t(scale(t(transposed_matrix)))
newnames <- lapply(
  rownames(scaled_transposed_matrix),
  function(x) bquote(italic(.(x))))
heatmp <- pheatmap(scaled_transposed_matrix, fontsize = 8,scale = "column",
                   fontsize_col = 8, cluster_rows = TRUE,
                   cluster_cols = FALSE, cellwidth = 10, cellheight = 10, angle_col = 45,
                   color = colorRampPalette(c("blue", "white", "red"))(50),
                   labels_row = as.expression(newnames))
png(filename = "Figures/Fig5/P21_ZI_heatmap_Dlx_DEG_Italic.png", width = 5000, height = 2000, res = 300)
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
#png(filename = "Figures/Fig5/P21_ZI_heatmap_Dlx_DEG_Italic.png", width = 3000, height = 2000)
#print(p)
#dev.off()

#p <- pheatmap(Avg.Scaling$SCT[valid_genes,], fontsize = 8, scale = "row",
#              fontsize_col = 8, cluster_rows = T,
#              cluster_cols = F, cellwidth = 10, cellheight = 10, angle_col = 45,
#              color = colorRampPalette(c("blue", "white", "red"))(50))
#png(filename = "Figures/Fig5/P21_ZI_heatmap_Dlx_DEG.png", width = 3000, height = 2000)
#print(p)
#dev.off()

#distribution
library(viridis)
library(viridisLite)
library(ggplot2)
Idents(DLXCKO) <- "Cluster"
identities <- Idents(DLXCKO)

#Grid <- DLXCKO@reductions$umap@cell.embeddings
#Grid <- data.frame(Grid)
#Grid$orig.ident <- DLXCKO$Genotype_Final
#Grid_ds <- Grid[,1:2]
#p <- ggplot(Grid, aes(x = umap_1, y = umap_2)) +
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
#png(filename = "Figures/Fig5/P21_ZI_densityheatmap_Dlx.png", width = 1200, height = 1000)
#print(p)
#dev.off()

#distribution 2
library(cetcolor)
library(Seurat)
library(ggplot2)
# generate UMAP plot
#thresholding
# Generate UMAP plot
pl1 <- DimPlot(DLXCKO, combine = F, cols = group, split.by = "Genotype_Final")
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
  stat_density_2d(aes_string(x = "umap_1", y = "umap_2", fill = "threshold_density(after_stat(level))"), 
                  linewidth = 0.2, geom = "density_2d_filled", 
                  colour = "ivory", alpha = 0.4, n = 150, h = c(1.2, 1.2)) & 
  scale_fill_gradientn(colours = scale.col, limits = c(0, 100)) & 
  DarkTheme()
ggsave(filename = "Figures/Fig5/P21_ZI_densityheatmap_Dlx2.eps",
       plot = gplot, width = 20, height = 10, units = "in", device = cairo_ps)
#png(filename = "Figures/Fig5/P21_ZI_densityheatmap_Dlx2.png", width = 1200, height = 1000)
#print(pl2)
#dev.off()

#density
library(dplyr)
library(reshape2)
library(Rmisc)
Idents(DLXCKO) <- "Cluster"
identities <- Idents(DLXCKO)
new_order <- sort(levels(identities))
Idents(DLXCKO) <- factor(identities, levels = new_order)
meta.data <- DLXCKO@meta.data
graph <- 100*prop.table(table(Idents(DLXCKO), DLXCKO$Genotype_Final), margin = 1) #margin 1= row, 2 = column
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
ggsave(filename = "Figures/Fig5/P21_ZI_densitypercentage_Dlx.eps",
       plot = gplot, width = 20, height = 10, units = "in", device = cairo_ps)
#png(filename = "Figures/Fig5/P21_ZI_densitypercentage_Dlx.png", width = 1200, height = 1000)
#print(p)
#dev.off()



####identify cluster####
Cells1 <-WhichCells(DLXCKO, idents = "AntID_ID_TT")
DimPlot(DLXCKO, reduction = "umap", label = F, 
        pt.size = 0.5, cols = , cells.highlight = list(Cells1)) +
  ggplot2::scale_color_manual(labels = c("Rest","Cluster"), values = c("lightgrey","#b54ccf"))


