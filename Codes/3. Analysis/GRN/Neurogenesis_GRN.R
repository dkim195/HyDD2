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

####Neurogenesis processing####
data <- read.table("SCENIC/eRegulon/neurogenesis_eRegulons_extended.tsv", header = TRUE, sep = "\t")
#data <- read.table("SCENIC/eRegulon/pattern_neurogenesis_eRegulons_extended.tsv", header = TRUE, sep = "\t")

# Filter and prepare the data
filtered_data <- data %>%
  filter(grepl("extended_\\+/\\+", eRegulon_name) | grepl("extended_\\-/\\+", eRegulon_name)) %>%
  select(eRegulon_name, Gene, importance_TF2G, triplet_rank) %>%
  mutate(type = ifelse(grepl("\\+/\\+", eRegulon_name), "Activator", ifelse(grepl("\\-/\\+", eRegulon_name), "Repressor", NA)))

selected_genes <- c("Agrp","Npy","Pomc","Gal",
                    "Pnoc","Th","Penk","Ghrh","Kiss1",
                    "Hcrt","Pmch","Sst","Pvalb","Tac1","Tac2",
                    "Grp","Avp","Oxt", "Trh", "Crh", "Vip", "Hdc",
                    "Slc32a1", "Slc17a6", "Slc17a7", "Cartpt", "Cck", "Nts") #Tub

selected_genes <- unique(selected_genes)
filtered_data <- filtered_data %>%
  filter(Gene %in% selected_genes)

filtered_data <- filtered_data %>%
  mutate(eRegulon_name = str_replace(eRegulon_name, "_extended_\\+/\\+", ""),
         eRegulon_name = str_replace(eRegulon_name, "_extended_\\-/\\+", ""))

selected_eRegulons <- c("Irx6","Irx5","Irx3","Lmx1b","Foxa1","Lmx1a","Foxa2","Barhl1",
                        "Bhlhe23", "Neurod6", "Ebf2", "Neurod4", "Nr4a2", "Nhlh1", "Neurod2", "Ebf1","Barhl2",
                        "Foxb1", "Pcp4", "Lhx1", "Fezf2", "Sim1", "Nhlh2", "Lhx5", "Neurod1", "Neurod2", "Uncx", "Nhlh1", "Emx2", "Nkx6-2",
                        "Otp", "Bsx", "Sim1", "Arxes2", "Bhlhe22", "NR5a1", "Sox14", "Nr0b1", "Satb2", "Fezf1",
                        "Neurog3","Six6","Rax","Sox3", "Hmx2", "Hmx3", "Gsx1", "Prox1", "Bsx", "Cited1",
                        "Meis2", "Sp8", "Sp9", "Arx", "Pax6", "Dlx6","Dlx1","Dlx2","Dlx5","Lhx6", "Onecut3",
                        "Isl1","Lhx2","Nkx2-1","Nkx2-2")

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
#p <- ggraph(ggraph_layout) + 
#  geom_edge_link(aes(color = type), 
#                 arrow = arrow(length = unit(filtered_data$normalized_rank, 'mm'), type = "closed"), 
#                 end_cap = circle(3, 'mm')) +
#  geom_point(aes(x = x, y = y), size = -1) +  # Use node positions for points
#  geom_text_repel(aes(x = x, y = y, label = name), size = 5, max.overlaps = 200) +  # Use geom_text_repel with node positions and increased max.overlaps
#  scale_edge_color_manual(values = c("Activator" = "blue", "Repressor" = "red")) +
#  scale_edge_width(range = c(0.5, 2)) +
#  theme_void() +
#  labs(title = "Pathway Plot of Selected eRegulons",
#       edge_color = "Type",
#       edge_width = "Importance TF2G",
#       edge_alpha = "Importance TF2G")

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

postscript(file = "Figures/Fig1/Neurogenesis_SCENIC_GRN/Neurogenesis_SCENIC_GRN.eps", width = 30, height = 15,
           horizontal = FALSE, onefile = FALSE, paper = "special")
print(p)
dev.off()

####All####
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

####Neurogenesis processing####
data <- read.table("SCENIC/eRegulon/neurogenesis_eRegulons_extended.tsv", header = TRUE, sep = "\t")
data <- read.table("SCENIC/eRegulon/pattern_neurogenesis_eRegulons_extended.tsv", header = TRUE, sep = "\t")
# Filter and prepare the data
filtered_data <- data %>%
  filter(grepl("extended_\\+/\\+", eRegulon_name) | grepl("extended_\\-/\\+", eRegulon_name)) %>%
  select(eRegulon_name, Gene, importance_TF2G, triplet_rank) %>%
  mutate(type = ifelse(grepl("\\+/\\+", eRegulon_name), "Activator", ifelse(grepl("\\-/\\+", eRegulon_name), "Repressor", NA)))

selected_genes <- c("Agrp","Npy","Pomc","Gal",
                    "Pnoc","Th","Penk","Ghrh","Kiss1",
                    "Hcrt","Pmch","Sst","Pvalb","Tac1","Tac2",
                    "Grp","Avp","Oxt", "Trh", "Crh", "Vip", "Hdc",
                    "Slc32a1", "Slc17a6", "Slc17a7", "Cartpt", "Cck", "Nts") #Tub

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
#p <- ggraph(ggraph_layout) + 
#  geom_edge_link(aes(color = type), 
#                 arrow = arrow(length = unit(filtered_data$normalized_rank, 'mm'), type = "closed"), 
#                 end_cap = circle(3, 'mm')) +
#  geom_point(aes(x = x, y = y), size = -1) +  # Use node positions for points
#  geom_text_repel(aes(x = x, y = y, label = name), size = 5, max.overlaps = 200) +  # Use geom_text_repel with node positions and increased max.overlaps
#  scale_edge_color_manual(values = c("Activator" = "blue", "Repressor" = "red")) +
#  scale_edge_width(range = c(0.5, 2)) +
#  theme_void() +
#  labs(title = "Pathway Plot of Selected eRegulons",
#       edge_color = "Type",
#       edge_width = "Importance TF2G",
#       edge_alpha = "Importance TF2G")

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


postscript(file = "Figures/Fig1/Neurogenesis_SCENIC_GRN/Neurogenesis_SCENIC_GRN2.eps", width = 30, height = 15,
           horizontal = FALSE, onefile = FALSE, paper = "special")
print(p)
dev.off()
