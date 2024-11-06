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

####neurogenesis####
load(file = "scRNA/Robj/neurogenesis.Robj")

####neurogenesis-scRNA clusters####
Idents(neurogenesis) <- "Cluster_Pass2"
neurogenesis <- subset(neurogenesis, idents = c("ARC (Npy_Agrp)", "ARC (Pomc)", "ARC_VMH", "ARC_VMH (Tac2_PMN)", 
                                               "DMH-LH (Grp_Cck)", "DMH-LH (Npvf)", "DMH (Npw_PMN)", "Isl1 cells", 
                                                "LH (Hcrt_Oxt)", "LH (Pmch_Trh)", "MMN", "MMN (Cck)", "MMN (Npy_Cck)", 
                                                "MMN (Nts_Tac1_Cck)", "PMN (Ghrh_Gal_Cited1)", "PMN (Ghrh_Gal_Oxt)", "POA_SCN", 
                                                "PreThal", "PreThal_ID", "PVN_SON (Avp_Gal_Oxt)", "SCN (Rorb)", 
                                                "SMN", "SMN (Calb2)", "SMN (Tac2)", "Sst_Npy", 
                                                "TMN_PH (Hdc)"))
neurogenesis <- SCTransform(neurogenesis,
                       vars.to.regress = c("nCount_RNA","nFeature_RNA"))
neurogenesis <- RunPCA(neurogenesis, npcs = 50, ndims.print = NA, verbose = F)
neurogenesis <- RunHarmony(neurogenesis, group.by.vars = "orig.ident",
                      assay.use = "SCT", plot.convergence = T)
neurogenesis <- RunUMAP(neurogenesis, reduction = "harmony", dims = 1:50,
                   n.neighbours = 10L, min.dist = 0.01, spread = 1)
DimPlot(neurogenesis, reduction = "umap", label = T,pt.size = 0.1) + NoLegend() + NoAxes()
save(neurogenesis, file = "scRNA/Robj/neurogenesis_neuron.Robj")


genes <- c("Hcrt","Kiss1","Pmch","Sst","Crh") # not detected in SCENIC, HCRT Rfx4 but 

genes <- c("Agrp", "Tac2", "Vip", "Avp", "Oxt", "Npy", "Gal",
           "Pomc","Th","Ghrh","Nts", "Pnoc", "Cartpt","Cck","Trh","Tac1",
           "Grp","Penk","Hdc") #SCENIC detection

genes <- c("Agrp","Npy","Gal","Pomc","Th","Ghrh","Nts",
           "Pnoc", "Cartpt", "Cck", "Trh", "Tac1", "Penk", "Hdc") #SCENIC & Xenium & ARC/TMN/MM

VlnPlot(neurogenesis, genes, stack = T, flip = T) + NoLegend()

####PLOT####
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
selected_genes <- c("Agrp","Npy","Gal","Pomc","Th","Ghrh","Nts",
                    "Pnoc", "Cartpt", "Cck", "Trh", "Tac1", "Penk", "Hdc") #SCENIC & Xenium & ARC/TMN/MM

filtered_data <- data %>%
  filter(grepl("extended_\\+/\\+", eRegulon_name) | grepl("extended_\\-/\\+", eRegulon_name)) %>%
  select(eRegulon_name, Gene, importance_TF2G, triplet_rank) %>%
  mutate(type = ifelse(grepl("\\+/\\+", eRegulon_name), "Activator", ifelse(grepl("\\-/\\+", eRegulon_name), "Repressor", NA)))

selected_genes <- unique(selected_genes)
filtered_data <- filtered_data %>%
  filter(Gene %in% selected_genes)

filtered_data <- filtered_data %>%
  mutate(eRegulon_name = str_replace(eRegulon_name, "_extended_\\+/\\+", ""),
         eRegulon_name = str_replace(eRegulon_name, "_extended_\\-/\\+", ""))

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
p

postscript(file = "Figures/Fig1/Neurogenesis_SCENIC_GRN/Neurogenesis_SCENIC_GRN_Neuropeptide_Xenium.eps", width = 30, height = 15,
           horizontal = FALSE, onefile = FALSE, paper = "special")
print(p)
dev.off()


####Plot2####
selected_genes <- c("Npy","Gal","Nts",
                    "Pnoc", "Cartpt", "Cck", "Tac1") #MM/smn

selected_eRegulons <- c("Foxb1", "Pcp4", "Lhx1", "Sim1", "Lhx5", "Nr4a2", "Lmx1a", "Lmx1b", "Barhl1", "Barhl2", "Ebf2",
                        "Nhlh1", "Nhlh2", "Ebf3", "Nkx2-4","Uncx")

filtered_data <- data %>%
  filter(grepl("extended_\\+/\\+", eRegulon_name) | grepl("extended_\\-/\\+", eRegulon_name)) %>%
  select(eRegulon_name, Gene, importance_TF2G, triplet_rank) %>%
  mutate(type = ifelse(grepl("\\+/\\+", eRegulon_name), "Activator", ifelse(grepl("\\-/\\+", eRegulon_name), "Repressor", NA)))

selected_genes <- unique(selected_genes)
filtered_data <- filtered_data %>%
  filter(Gene %in% selected_genes)

filtered_data <- filtered_data %>%
  mutate(eRegulon_name = str_replace(eRegulon_name, "_extended_\\+/\\+", ""),
         eRegulon_name = str_replace(eRegulon_name, "_extended_\\-/\\+", ""))

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
p

postscript(file = "Figures/Fig1/Neurogenesis_SCENIC_GRN/Neurogenesis_SCENIC_GRN_Neuropeptide_Xenium_MM.eps", width = 30, height = 15,
           horizontal = FALSE, onefile = FALSE, paper = "special")
print(p)
dev.off()