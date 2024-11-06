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
library(grid)
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram
setwd("/media/thomaskim/Data/")

####pattern####
load(file = "scRNA/Robj/pattern.Robj")

####pattern-scRNA clusters####
#POA/SC -> PVH/SON
#PVH/SON ANTVEN = PVH/SON
#PVH/SON = PVH/SON (PRO)
#SST(DMH) = PRETHAL (SP8)
#PRETHAL/ADNID = PRETHAL (SP8)
#ALL TUB LH, ARC/VMH, PMCH = TUB
Idents(pattern) <- "Cluster_Pass2"
pattern <- RenameIdents(pattern, "Ant_ID" = "AntID_ID",
                        "NPC (gliogenic?)" = "NPC (gliogenic)",
                        "PMN? (Gsx1)" = "PMN (Gsx1)",
                        "PVH_SON" = "PVH_SON (Pro)",
                        "PreThal_AntID" = "PreThal (Sp8)",
                        "SST (DMH?)" = "PreThal (Sp8)",
                        "POA/SCN?" = "PVH_SON",
                        "PVH_SON_AntVen" = "PVH_SON",
                        "ZLI?" = "ZLI",
                        "Tub (LH-Pmch)" = "Tub",
                        "Tub (LH)" = "Tub",
                        "Tub (ARC_VMH)" = "Tub")
identities <- Idents(pattern)
new_order <- sort(levels(identities))
Idents(pattern) <- factor(identities, levels = new_order)
group <- c(
  "AntID_ID" = "#808AA9", 
  "MMN" = "#F0A881", 
  "MMN (Prog)" = "#F1CCB8",
  "Neural Pro (Ascl1)" = "#00FF1B",  
  "Neural Pro (Neurog2)" = "#4EF860",
  "NPC" = "#F6B4F1", 
  "NPC (Ascl1)" = "#95D992",
  "NPC (gliogenic)" = "#E6B3B3", 
  "NPC (Neurog2)" = "#C1E6BF",
  "NSC" = "#AC8DC4", 
  "NSC (Neurog2)" = "#C3AED3",
  "NSC_NPC" = "#4EFAA8", 
  "PMN" = "#FFBB40", 
  "PMN (Gsx1)" = "#DD8F00", "PMN (Pro)" = "#DC7C7C",
  "PreThal (Pro)" =  "#D3D3D3",
  "PreThal (Sp8)" = "#3A62D9",
  "PreThal (Sp9)" = "#979797",
  "PVH_SON (Pro)" = "#DFCE93", 
  "PVH_SON" = "#E9DEC1",
  "SMN" = "#93F6F8",
  "SMN (Prog)" = "#A7D4E7",
  "Tub"= "#FF4136",
  "Tub (Prog)" = "#CE6ACB", 
  "ZLI"="#434D57")

#UMAP
p <- DimPlot(pattern, label = T, cols = group, raster = F)
png(filename = "Figures/Fig1/Pattern_scRNA/DEG_pattern_scRNA_raw.png", width = 900, height = 900)
print(p)
dev.off()

#UMAP
p <- DimPlot(pattern, label = , cols = group, raster = F)+ NoLegend() + NoAxes()
png(filename = "Figures/Fig1/Pattern_scRNA/DEG_pattern_scRNA.png", width = 3000, height = 3000, res = 300)
print(p)
dev.off()

markers <- FindAllMarkers(pattern, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "Figures/Fig1/Pattern_scRNA/DEG_pattern_scRNA.csv")

#plot
library(pheatmap)
lists <- read.csv(file = "Figures/Fig1/Pattern_scRNA/DEG_pattern_scRNA.csv",
                  header = T, row.names = 1)
#selected groups
select_group <- c("Neural Pro (Neurog2)", "MMN (Prog)", "MMN", "SMN (Prog)","SMN",
                  "PVH_SON (Pro)", "PVH_SON", "Neural Pro (Ascl1)", "Tub (Prog)", "Tub",
                  "PMN (Pro)", "PMN", "PreThal (Pro)", "PMN (Gsx1)", "PreThal (Sp9)", "AntID_ID", "PreThal (Sp8)")
# Load transcription factor list
library(data.table)
tf_list_url <- "http://humantfs2.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt"
tf_list <- fread(tf_list_url, header = FALSE, sep = "\t")
transcription_factors <- tf_list$V1
# Function to capitalize only the first character
capitalize_first <- function(gene) {
  paste0(toupper(substr(gene, 1, 1)), tolower(substr(gene, 2, nchar(gene))))
}
# Apply the function to the list of genes
formatted_genes <- sapply(transcription_factors, capitalize_first)
lists <- lists[lists$gene %in% formatted_genes, ]

genes <- lists %>%
  filter(cluster %in% select_group) %>%  # Keep only select_group
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(avg_log2FC > 3) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  distinct(gene) %>%  # Ensure unique genes are selected if they appear in multiple clusters
  pull(gene)  # Extract the gene names as a character vector

#genes <- lists %>% arrange(-avg_log2FC) %>% filter(avg_log2FC> 3)
#genes<- as.character(genes$gene)
Scaling  <- ScaleData(pattern, features = genes)
Avg.Scaling <- AverageExpression(Scaling, assay = "SCT", slot = "counts",
                                 verbose = T) #exponential minus 1 mean(expm1(x)) of the SCTTransform value
transposed_matrix <- t(Avg.Scaling$SCT[genes,])
scaled_transposed_matrix <- t(scale(t(transposed_matrix)))
newnames <- lapply(
  rownames(scaled_transposed_matrix),
  function(x) bquote(italic(.(x))))
heatmp <- pheatmap(scaled_transposed_matrix, fontsize = 8,scale = "column",
                   fontsize_col = 8, cluster_rows =F,
                   cluster_cols = F, cellwidth = 10, cellheight = 10, angle_col = 45,
                   color = colorRampPalette(c("blue", "white", "red"))(50),
                   labels_row = as.expression(newnames))
library(grDevices)
eps_file <- "Figures/Fig1/Pattern_scRNA/Heatmap_pattern_scRNA_Italic.eps"
setEPS()
postscript(eps_file, width = 8, height = 15) 
print(heatmp)
dev.off()


identities <- Idents(pattern)
new_order <- sort(levels(identities))
Idents(pattern) <- factor(identities, levels = new_order)
meta.data <- pattern@meta.data
graph <- 100*prop.table(table(Idents(pattern), pattern$Age_Sum), margin = 1) #margin 1= row, 2 = column
dat <- melt(graph)
write.csv(dat, file = "Figures/Fig1/Pattern_scRNA/dat.csv", row.names = F)
#distribution
library(viridis)
library(viridisLite)
library(ggplot2)
library(colorspace)
library(RColorBrewer)
library(cetcolor)
pl1 <- DimPlot(pattern, combine = F, cols = group, split.by = "Age_Sum", raster = F, pt.size = 0.01)
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
  DarkTheme() & 
  theme(legend.position = "none") # Remove legend
ggsave(filename = "Figures/Fig1/Pattern_scRNA/densityheatmap_pattern.eps",
       plot = gplot, width = 20, height = 10, units = "in", device = cairo_ps)




####pattern-matching scATAC####
Idents(pattern) <- "Cluster_Pass2"
pattern <- RenameIdents(pattern, "Ant_ID" = "AntID_ID", 
                        "MMN (Prog)" = "MMN", 
                        "PMN (Pro)" = "PMN",
                        "PMN? (Gsx1)" = "PMN",
                        "PreThal (Pro)" = "PreThal", 
                        "PreThal (Sp8)" = "PreThal",
                        "PreThal (Sp9)" = "PreThal",
                        "PreThal_AntID" = "PreThal",
                        "PVH_SON_AntVen" = "PVH_SON",
                        "SMN (Prog)" = "SMN",
                        "Tub (ARC_VMH)" = "Tub", 
                        "Tub (LH-Pmch)"= "Tub",
                        "Tub (LH)"= "Tub",
                        "Tub (Prog)"= "Tub",
                        "NPC (gliogenic?)" = "NPC (gliogenic)",
                        "PMN? (Gsx1)" = "PMN (Gsx1)",
                        "POA/SCN?"= "PVH_SON",
                        "SST (DMH?)" = "PreThal", 
                        "ZLI?" = "ZLI")

# Reorder levels alphabetically
identities <- Idents(pattern)
new_order <- sort(levels(identities))
Idents(pattern) <- factor(identities, levels = new_order)
group <- c("AntID_ID" = "#808AA9", 
           "MMN"= "#F0A881",  "Neural Pro (Ascl1)" = "#00FF1B", 
           "Neural Pro (Neurog2)" = "#4EF860", "NPC" = "#F6B4F1", "NPC (Ascl1)" = "#95D992", 
           "NPC (gliogenic)" = "#E6B3B3", "NPC (Neurog2)" = "#C3AED3", "NSC" = "#AC8DC4",
           "NSC (Neurog2)" = "#C3AED3", "NSC_NPC" = "#4EFAA8", 
           "PMN" = "#FFBB40", "PreThal" = "#3A62D9",
           "PVH_SON"= "#DFCE93", "SMN"= "#93F6F8", "Tub"= "#FF4136", 
           "ZLI"="#434D57")

#UMAP
p <- DimPlot(pattern, label = T, cols = group, raster = F)
png(filename = "Figures/Fig1/Pattern_scRNA/UMAP_pattern_scRNA_compact_raw.png", width = 900, height = 900)
print(p)
dev.off()

#UMAP
p <- DimPlot(pattern, label = , cols = group, raster = F)+ NoLegend() + NoAxes()
png(filename = "Figures/Fig1/Pattern_scRNA/UMAP_pattern_scRNA_compact.png", width = 3000, height = 3000, res = 300)
print(p)
dev.off()

markers <- FindAllMarkers(pattern, test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "Figures/Fig1/Pattern_scRNA/DEG_pattern_scRNA_compact.csv")

library(pheatmap)
lists <- read.csv(file = "Figures/Fig1/Pattern_scRNA/DEG_pattern_scRNA_compact.csv",
                  header = T, row.names = 1)
# Load transcription factor list
library(data.table)
tf_list_url <- "http://humantfs2.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt"
tf_list <- fread(tf_list_url, header = FALSE, sep = "\t")
transcription_factors <- tf_list$V1
# Function to capitalize only the first character
capitalize_first <- function(gene) {
  paste0(toupper(substr(gene, 1, 1)), tolower(substr(gene, 2, nchar(gene))))
}
# Apply the function to the list of genes
formatted_genes <- sapply(transcription_factors, capitalize_first)
lists <- lists[lists$gene %in% formatted_genes, ]

#all groups
genes <- lists %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(avg_log2FC > 3) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  distinct(gene) %>%  # Ensure unique genes are selected if they appear in multiple clusters
  pull(gene)  # Extract the gene names as a character vector

Scaling  <- ScaleData(pattern, features = genes)
Avg.Scaling <- AverageExpression(Scaling, assay = "SCT", slot = "counts",
                                 verbose = T) #exponential minus 1 mean(expm1(x)) of the SCTTransform value
transposed_matrix <- t(Avg.Scaling$SCT[genes,])
scaled_transposed_matrix <- t(scale(t(transposed_matrix)))
newnames <- lapply(
  rownames(scaled_transposed_matrix),
  function(x) bquote(italic(.(x))))
heatmp <- pheatmap(scaled_transposed_matrix, fontsize = 8,scale = "column",
                   fontsize_col = 8, cluster_rows = TRUE,
                   cluster_cols = FALSE, cellwidth = 10, cellheight = 10, angle_col = 45,
                   color = colorRampPalette(c("blue", "white", "red"))(50),
                   labels_row = as.expression(newnames))
library(grDevices)
eps_file <- "Figures/Fig1/Pattern_scRNA/Heatmap_pattern_scRNA_compact_Italic.eps"
setEPS()
postscript(eps_file, width = 8, height = 15) 
print(heatmp)
dev.off()

