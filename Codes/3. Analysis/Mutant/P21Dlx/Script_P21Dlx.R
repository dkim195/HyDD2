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
load(file = "scRNA/Robj/DLXCKO_P21_ZITRN_Neuron_Final.Robj")

identities <- Idents(DLXCKO)
new_order <- sort(levels(identities))
Idents(DLXCKO) <- factor(identities, levels = new_order)
meta.data <- DLXCKO@meta.data
graph <- 100*prop.table(table(Idents(DLXCKO), DLXCKO$Genotype_Final), margin = 1) #margin 1= row, 2 = column
dat <- melt(graph)
write.csv(dat, file = "Figures/P21Dlx/dat.csv", row.names = F)


#Find unique regions
DLXCKO <- FindNeighbors(DLXCKO, dims = 1:20, reduction = "harmony")
DLXCKO <- FindClusters(DLXCKO, resolution = 0.8) #SCT_snn_res.0.8
DimPlot(DLXCKO, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

contingency_table <- table(DLXCKO@active.ident, DLXCKO$Genotype_Final)
percentage_table <- prop.table(contingency_table, margin = 1) * 100

library(tidyr)
library(dplyr)
reshaped_df <- as.data.frame(percentage_table) %>%
  pivot_wider(names_from = Var2, values_from = Freq, id_cols = Var1)
reshaped_df$Mut_Control_Ratio <- with(reshaped_df, Mut / Ctrl)
View(reshaped_df)

clusters <-c("14","17","11","20","2","19","13","18","3","16","6","1","15","8")
for (cluster in clusters) {
  # Find markers for the current cluster
  markers <- FindMarkers(DLXCKO, ident.1 = cluster,  # No need to convert if already a character
                         logfc.threshold = 0.5,
                         min.pct = 0.2, verbose = TRUE)
  
  # Define the file path using %s for string formatting
  file_name <- sprintf("Figures/P21Dlx/P21Dlx_Enriched_DLXCKO_C%s.csv", cluster)
  
  # Save to CSV
  write.csv(markers, file_name, row.names = T)
}

#p <- FeaturePlot_scCustom(seurat_object = DLXCKO, features = c("Lhx9"),
#                          alpha_exp = 0.75, split.by = "Genotype") & NoLegend() & NoAxes() +
#  theme(plot.title = element_blank())
#png(filename = "Figures/Dlx/Lhx9.png", width = 3000, height = 3000, res = 300)
#print(p)
#dev.off()

library(scCustomize)
plot_and_save_genes <- function(seurat_object, features, directory = "Figures") {
  # Ensure the directory exists
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  for (gene in features) {
    # Generate the plot for the current gene
    p <- FeaturePlot_scCustom(seurat_object = seurat_object, features = c(gene),
                              order = T, split.by = "Genotype_Final", pt.size = 0.1, na_cutoff = 1) & NoLegend() & NoAxes() +
      theme(plot.title = element_blank())
    
    # Define the file path
    file_path <- file.path(directory, paste0(gene, ".png"))
    
    # Save the plot to a PNG file
    png(filename = file_path, width = 2400, height = 1200, res = 300)
    print(p)
    dev.off()
  }
}

genes <- c("Oxt","Esr2","Sim1","Reln","Otp","Hcrt","Sst","Bsg",
           "Pnoc","Penk","Pvalb","Cartpt","Pmch","Atoh7","Gal","Tbr1","Crh","Avp","Trh","Lef1","Tac1","Spp1","Ecel1")
Stacked_VlnPlot(seurat_object = DLXCKO, features = genes, x_lab_rotate = TRUE, group.by = "Cluster")

VlnPlot(DLXCKO, genes, group.by = "Cluster")
plot_and_save_genes(DLXCKO, genes, "Figures/P21Dlx")


####GRN####
library(data.table)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(ggrepel)
library(reshape2)
library(pheatmap)
library(stringr)
# Set working directory and load the data
setwd("/media/thomaskim/Data/")
data <- read.table("SCENIC/eRegulon/pattern_neurogenesis_eRegulons_extended.tsv", header = TRUE, sep = "\t")

# Filter and prepare the data
filtered_data <- data %>%
  filter(grepl("extended_\\+/\\+", eRegulon_name) | grepl("extended_\\-/\\+", eRegulon_name)) %>%
  select(eRegulon_name, Gene, importance_TF2G, triplet_rank) %>%
  mutate(type = ifelse(grepl("\\+/\\+", eRegulon_name), "Activator", ifelse(grepl("\\-/\\+", eRegulon_name), "Repressor", NA)))

selected_genes <- c("Meis2","Arxes1","Arxes2","Spp1","Ecel1","Tac1","Pvalb","Penk","Sst","Gal") #Tub

selected_genes <- unique(selected_genes)
filtered_data <- filtered_data %>%
  filter(Gene %in% selected_genes)

filtered_data <- filtered_data %>%
  mutate(eRegulon_name = str_replace(eRegulon_name, "_extended_\\+/\\+", ""),
         eRegulon_name = str_replace(eRegulon_name, "_extended_\\-/\\+", ""))

selected_eRegulons <- c("Dlx1","Dlx2")

filtered_data <- filtered_data %>%
  filter(eRegulon_name %in% selected_eRegulons)

# Convert eRegulon_name to a factor with levels ordered by the original eRegulon_name
filtered_data <- filtered_data %>%
  mutate(eRegulon_name = factor(eRegulon_name, levels = unique(eRegulon_name)))

filtered_data <- filtered_data %>%
  mutate(normalized_rank = (max(triplet_rank) - triplet_rank) / max(triplet_rank) * 4 + 1)  # Scale to range 1-5 for arrow length

graph_gene <- graph_from_data_frame(filtered_data %>% select(eRegulon_name,
                                                             Gene, type, 
                                                             importance_TF2G, normalized_rank), directed = TRUE)
ggraph_layout <- create_layout(graph_gene, layout = 'fr')

p <- ggraph(ggraph_layout) + 
  geom_edge_link(aes(color = type), 
                 arrow = arrow(length = unit(3, 'mm'), type = "closed"),  # Set a constant length for all arrows
                 end_cap = circle(3, 'mm'),
                 edge_width = 1) +  # Set a constant width for all edges
  geom_point(aes(x = x, y = y), size = -1) +  # Use node positions for points
  geom_text_repel(aes(x = x, y = y, label = name), size = 5, max.overlaps = 200) +  # Use geom_text_repel with node positions and increased max.overlaps
  scale_edge_color_manual(values = c("Activator" = "blue", "Repressor" = "red")) +
  theme_void() +
  labs(title = "Pathway Plot of Selected eRegulons (Top 10 Percentile by Triplet Rank)",
       edge_color = "Type")

print(p)
postscript(file = "Figures/P21Dlx/SCENIC_Dlx_eRegulon_Pathway_Regulon_TF.eps", width = 30, height = 15,
           horizontal = FALSE, onefile = FALSE, paper = "special")
print(p)
dev.off()