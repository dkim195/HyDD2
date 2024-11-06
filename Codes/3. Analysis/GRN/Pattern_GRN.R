library(data.table)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(ggrepel)
library(reshape2)
library(pheatmap)
library(stringr)

# Set working directory
setwd("/media/thomaskim/Data/")

# Load transcription factor list
tf_list_url <- "http://humantfs2.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt"
tf_list <- fread(tf_list_url, header = FALSE, sep = "\t")
transcription_factors <- tf_list$V1

# Function to capitalize only the first character
capitalize_first <- function(gene) {
  paste0(toupper(substr(gene, 1, 1)), tolower(substr(gene, 2, nchar(gene))))
}

# Apply the function to the list of genes
formatted_genes <- sapply(transcription_factors, capitalize_first)

# Load eRegulon data
data <- read.table("SCENIC/eRegulon/PATTERN_eRegulons_extended.tsv", header = TRUE, sep = "\t")

# Filter and prepare the data
filtered_data <- data %>%
  filter(grepl("extended_\\+/\\+", eRegulon_name) | grepl("extended_\\-/\\+", eRegulon_name)) %>%
  select(eRegulon_name, Gene, importance_TF2G, triplet_rank) %>%
  mutate(type = ifelse(grepl("\\+/\\+", eRegulon_name), "Activator", ifelse(grepl("\\-/\\+", eRegulon_name), "Repressor", NA)))

# Select specific eRegulons
selected_eRegulons <- c("Isl1_extended_+/+", "Isl1_extended_-/+")

filtered_data <- filtered_data %>%
  filter(eRegulon_name %in% selected_eRegulons)

# Combine 'Isl1_extended_+/+' and 'Isl1_extended_-/+' into 'Isl1'
filtered_data <- filtered_data %>%
  mutate(eRegulon_name = str_replace(eRegulon_name, "_extended_\\+/\\+", ""),
         eRegulon_name = str_replace(eRegulon_name, "_extended_\\-/\\+", ""))

# Convert eRegulon_name to a factor with levels ordered by the original eRegulon_name
filtered_data <- filtered_data %>%
  mutate(eRegulon_name = factor(eRegulon_name, levels = unique(eRegulon_name)))

# Keep only transcription factors from the formatted list
filtered_data <- filtered_data[filtered_data$Gene %in% formatted_genes, ]

# Normalize the triplet_rank for arrow length mapping
filtered_data <- filtered_data %>%
  mutate(normalized_rank = (max(triplet_rank) - triplet_rank) / max(triplet_rank) * 4 + 1)  # Scale to range 1-5 for arrow length

# Load expression data to extract brain regions
lists <- read.csv(file = "Figures/Original/Fig1/Pattern_scRNA/DEG_pattern_scRNA_compact.csv",
                  header = T, row.names = 1)


desired_clusters <- c("AntID_ID", 
                      "MMN", "Neural Pro (Ascl1)", "Neural Pro (Neurog2)", "NPC", "NPC (Ascl1)", 
                      "NPC (gliogenic)", "NPC (Neurog2)", "NSC",
                      "PMN", "POA/SCN", "PreThal", "PVH_SON", "SMN", "Tub")

# Filter genes with avg_log2FC > 2
filtered_genes <- lists %>%
  filter(cluster %in% desired_clusters) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(avg_log2FC > 1)

# Extract the gene and corresponding brain region (cluster) information
gene_brain_regions <- filtered_genes %>%
  select(Gene = gene, BrainRegion = cluster)

# Merge brain region information into the filtered data
filtered_data <- merge(filtered_data, gene_brain_regions, by = "Gene", all.x = TRUE)

# Replace NA values in the BrainRegion column with "UNKNOWN"
filtered_data$BrainRegion[is.na(filtered_data$BrainRegion)] <- "UNKNOWN"

# Create graph object including BrainRegion as an attribute
graph_eRegulon <- graph_from_data_frame(filtered_data %>% 
                                          select(eRegulon_name, Gene, type, importance_TF2G, normalized_rank, BrainRegion), 
                                        directed = TRUE)

# Add BrainRegion attribute to graph vertices
V(graph_eRegulon)$BrainRegion <- filtered_data$BrainRegion[match(V(graph_eRegulon)$name, filtered_data$Gene)]

# Create the layout
ggraph_layout <- create_layout(graph_eRegulon, layout = 'fr')

# Plot with additional brain region separation
p <- ggraph(ggraph_layout) + 
  geom_edge_link(aes(color = type), 
                 arrow = arrow(length = unit(filtered_data$normalized_rank, 'mm'), type = "closed"), 
                 end_cap = circle(3, 'mm')) +
  geom_node_point(aes(x = x, y = y, color = BrainRegion), size = 3) +  # Color nodes by brain region
  geom_node_text(aes(x = x, y = y, label = name), repel = TRUE, size = 5, max.overlaps = 200) +  # Label nodes
  scale_edge_color_manual(values = c("Activator" = "blue", "Repressor" = "red")) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_color_manual(values = c("AntID_ID" = "#808AA9", 
                                "MMN"= "#F0A881",  "Neural Pro (Ascl1)" = "#00FF1B", 
                                "Neural Pro (Neurog2)" = "#4EF860", "NPC" = "#F6B4F1", "NPC (Ascl1)" = "#95D992", 
                                "NPC (gliogenic)" = "#E6B3B3", "NPC (Neurog2)" = "#C3AED3", "NSC" = "#AC8DC4",
                                "NSC (Neurog2)" = "#C3AED3", "NSC_NPC" = "#4EFAA8", 
                                "PMN" = "#FFBB40", "POA/SCN" = "#AA102E", "PreThal" = "#3A62D9",
                                "PVH_SON"= "#DFCE93", "SMN"= "#93F6F8", "SST (DMH)" = "#A700FF", "Tub"= "#FF4136", 
                                "ZLI"="#434D57", "UNKNOWN" = "gray")) +  # Define colors for regions
  theme_void() +
  labs(title = "Pathway Plot of Selected eRegulons (Top 10 Percentile by Triplet Rank) with Brain Regions",
       edge_color = "Type",
       edge_width = "Importance TF2G",
       color = "Brain Region")

# Print the plot
print(p)
