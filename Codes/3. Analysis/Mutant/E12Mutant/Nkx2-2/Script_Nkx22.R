####GRN PLot-NEW####
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

#only plot TF
#only plot Region-specific TF

# Load traG2M NPCription factor list
tf_list_url <- "http://humantfs2.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt"
tf_list <- fread(tf_list_url, header = FALSE, sep = "\t")
traG2M NPCription_factors <- tf_list$V1

# Function to capitalize only the first character
capitalize_first <- function(gene) {
  paste0(toupper(substr(gene, 1, 1)), tolower(substr(gene, 2, nchar(gene))))
}

# Apply the function to the list of genes
formatted_genes <- sapply(traG2M NPCription_factors, capitalize_first)

# Load eRegulon data
data <- read.table("SCENIC/eRegulon/PATTERN_eRegulons_extended.tsv", header = TRUE, sep = "\t")

# Filter and prepare the data
filtered_data <- data %>%
  filter(grepl("extended_\\+/\\+", eRegulon_name) | grepl("extended_\\-/\\+", eRegulon_name)) %>%
  select(eRegulon_name, Gene, importance_TF2G, triplet_rank) %>%
  mutate(type = ifelse(grepl("\\+/\\+", eRegulon_name), "Activator", ifelse(grepl("\\-/\\+", eRegulon_name), "Repressor", NA)))

# Select specific eRegulons
selected_eRegulons <- c("Nkx2-2_extended_+/+", "Nkx2-2_extended_-/+")

filtered_data <- filtered_data %>%
  filter(eRegulon_name %in% selected_eRegulons)

# Combine 'Nkx2-2_extended_+/+' and 'Nkx2-2_extended_-/+' into 'Nkx2-2'
filtered_data <- filtered_data %>%
  mutate(eRegulon_name = str_replace(eRegulon_name, "_extended_\\+/\\+", ""),
         eRegulon_name = str_replace(eRegulon_name, "_extended_\\-/\\+", ""))

# Convert eRegulon_name to a factor with levels ordered by the original eRegulon_name
filtered_data <- filtered_data %>%
  mutate(eRegulon_name = factor(eRegulon_name, levels = unique(eRegulon_name)))

# Keep only traG2M NPCription factors from the formatted list
filtered_data <- filtered_data[filtered_data$Gene %in% formatted_genes, ]

# Normalize the triplet_rank for arrow length mapping
filtered_data <- filtered_data %>%
  mutate(normalized_rank = (max(triplet_rank) - triplet_rank) / max(triplet_rank) * 4 + 1)  # Scale to range 1-5 for arrow length

# Load your CSV file
lists <- read.csv(file = "Figures/Original/Fig1/Pattern_scRNA/DEG_pattern_scRNA.csv",
                  header = TRUE, row.names = 1)
# Filter genes with avg_log2FC > 1 and group by Gene to identify shared genes
filtered_genes <- lists %>%
  filter(avg_log2FC > 1.2) %>%
  group_by(gene) %>%
  summarize(BrainRegion = paste(sort(unique(cluster)), collapse = "_")) %>%
  ungroup()
# If you want to view the resulting dataframe
filtered_genes <- filtered_genes[filtered_genes$gene %in% formatted_genes, ]
print(filtered_genes)
write.csv(filtered_genes, file = "Figures/Pattern_Gene.csv", row.names = F)

#Define
#1. if 1 gene = 2 Postmitotic Region, Region1_Region2 (i.e. Deaf1)
#2. if 1 gene = expressed from NPC/Neural Pro -> Single Region Progenitor, it is Region Progenitor marker (i.e. Bhlhe22)
#3. if 1 gene = expressed from NPC/Neural Pro -> Multiple Region Progenitor, it is NPC/Neural Pro marker (i.e Ascl1, Neurog2)
#4. if 1 gene = expressed from NPC/Neural Pro -> Region Progenitor -> Region(S), it is Regional marker and Rule#1 apply (i.e. Dlx1, Arid5b)
#5. some exceptions when the gene is shared by the regions that are thought to be independent (i.e. prethal and s/mammillary) and distant (pvh/son and s/mamm)
#6. G2M NPC excluded
#7. more than 1 Postmitotic region = Excluded
#8. PreThal Pro = PreThal as the development happens sooner than our collected data

MMN_Prog <- c("Akna","Cebpa","Emx2","Fezf2","Lhx5", "Nhlh2","Neurod1","Neurod2","Nkx6-2","Olig2","Uncx")
MMN <- c("Foxb1","Nkx2-4")
Tub <- c("Arid5b","Lhx9","Myc","Nr0b1","Nr5a1","Tead1","Satb2",
         "Lhx2","Klf6","Lcorl","Gbx2","Mef2c","Sox14","Sox5")
AntID_ID_PreThal <- c("Arx","Cux2","Dlx6","Etv1","Maf","Msantd1",
                      "Meis2","Pbx1","Zscan18","Zmat4","Zkscan2","Sp8","Sp9","Tbr1","Onecut3","Pax6","Ikzf4",
                      "Hivep2","Gsx2","Meis1","Lhx6","Nfil3","Csrnp3","Esrrg","Pou3f3","Vax1","Six3","Pou6f2","Nkx2-3","Zeb1","Zeb2","Tcf12")
NPC_Ascl1 <- c("Ascl1","Gli2","Nr2f6","Zic5","Terf1","Setbp1","Lin54","Nr2e1")
NPC_Neurog2 <-c("Atf3","Msx1","Tfap4","Smad3","Klf5","Klf2","Dbx1","Mecom","Neurog1","Neurog2","Nfe2l2",
                "Nfix","Otx1","Otx2","Plscr1","Prrx1","Sall4","Olig1","Rfx4")
PVH_SON <- c("Atf5","Cxxc5","Ddit3","Otp", "Xbp1","Insm2","Rxrg")
SMN <- c("Barhl1","Bcl11b","Dmrta2","Foxa1","Foxa2","Foxp1",
         "Foxp2","Irx1","Irx3","Irx5","Irx6","Lmx1a","Lmx1b","Myt1l","Pknox2","Pou3f4")
SMN_Prog <- c("Bhlhe23","Dach1","Tcf4","Ebf2","Ebf1","Ebf3","Nhlh1","Nr4a2","Barhl2","Neurod4")
Tub_PMN_Prog <- c("Fli1","Foxd2","Prdm13","Ptf1a","Zbtb20","Rax","Six6","Skor1","Neurog3")
PMN <- c("Mbnl2","Peg3","Hmx2","Hmx3","Lef1","Dach2", "Cited1", "Tcf7l2")
MMN_SMN <- c("Pitx2")
PVH_SON_Tub_PMN <- c("Fezf1","Bsx")
PVH_SON_MMN_SMN <- c("Dpf3","Neurod6","Sim1")
AntID_ID_PreThal_SMN <- c("Bcl11b","Myt1","Pbx3","Pou3f1","Zbtb38","Zfhx4","Scrt2","Tshz3")
AntID_ID_PreThal_PMN <- c("Zfpm2","Isl1","Gsx1","Mafb","Dlx1","Dlx2","Dlx5","Prox1","Rara","Tshz2")
AntID_ID_PreThal_MMN <- c("Zfhx3","Lhx1","Pou2f2")
AntID_ID_PreThal_PVH_SON <- c("Onecut1","Zbtb7c")

# Create a named list
gene_list <- list(
  SMN = SMN,
  SMN_Prog = SMN_Prog,
  MMN = MMN,
  MMN_Prog = MMN_Prog,
  Tub = Tub,
  PVH_SON = PVH_SON,
  PVH_SON_Tub_PMN = PVH_SON_Tub_PMN,
  PVH_SON_MMN_SMN = PVH_SON_MMN_SMN,
  Tub_PMN_Prog = Tub_PMN_Prog,
  PMN = PMN,
  AntID_ID_PreThal = AntID_ID_PreThal,
  AntID_ID_PreThal_SMN = AntID_ID_PreThal_SMN,
  AntID_ID_PreThal_PMN = AntID_ID_PreThal_PMN,
  AntID_ID_PreThal_MMN = AntID_ID_PreThal_MMN,
  AntID_ID_PreThal_PVH_SON = AntID_ID_PreThal_PVH_SON,
  MMN_SMN = MMN_SMN,
  NPC_Ascl1 = NPC_Ascl1,
  NPC_Neurog2 = NPC_Neurog2
)

# Transform the list into a dataframe
gene_brain_regions <- do.call(rbind, lapply(names(gene_list), function(region) {
  data.frame(Gene = gene_list[[region]], BrainRegion = region, stringsAsFactors = FALSE)
}))

# Merge brain region information into the filtered data
filtered_data <- merge(filtered_data, gene_brain_regions, by = "Gene", all.x = TRUE)

# Replace NA values in the BrainRegion column with "UNKNOWN"
filtered_data$BrainRegion[is.na(filtered_data$BrainRegion)] <- "UNKNOWN"

# Exclude UNKNOWN BrainRegion
filtered_data <- filtered_data %>% filter(BrainRegion != "UNKNOWN")

# Create graph object including BrainRegion as an attribute
graph_eRegulon <- graph_from_data_frame(filtered_data %>% 
                                          select(eRegulon_name, Gene, type, importance_TF2G, normalized_rank, BrainRegion), 
                                        directed = TRUE)

# Add BrainRegion attribute to graph vertices
V(graph_eRegulon)$BrainRegion <- filtered_data$BrainRegion[match(V(graph_eRegulon)$name, filtered_data$Gene)]

# Add type attribute to graph vertices
V(graph_eRegulon)$type <- filtered_data$type[match(V(graph_eRegulon)$name, filtered_data$Gene)]

# Create the layout
ggraph_layout <- create_layout(graph_eRegulon, layout = 'fr')

# Assign x and y positions based on type
for (i in seq_along(V(graph_eRegulon)$name)) {
  if (!is.na(V(graph_eRegulon)$type[i])) {  # Add a check for NA values
    if (V(graph_eRegulon)$type[i] == "Activator") {
      ggraph_layout$x[i] <- abs(ggraph_layout$x[i])  # Place activators on the right side
    } else if (V(graph_eRegulon)$type[i] == "Repressor") {
      ggraph_layout$x[i] <- -abs(ggraph_layout$x[i])  # Place repressors on the left side
    }
  }
}

# Get the unique brain regions present in the ggraph layout
present_brain_regions <- unique(V(graph_eRegulon)$BrainRegion)

# Define the full color palette for all possible brain regions
full_color_palette <- c(
  "MMN"= "#F0A881",  
  "SMN"= "#93F6F8", 
  "SMN_Prog" = "#4EF860", 
  "MMN_Prog" = "#95D992", 
  "Tub" = "#FF4136", 
  "PVH_SON" = "#DFCE93", 
  "PVH_SON_Tub_PMN" = "#AA102E", 
  "PVH_SON_MMN_SMN" = "#A700FF", 
  "Tub_PMN_Prog" = "#C3AED3", 
  "PMN" = "#FFBB40", 
  "AntID_ID_PreThal" = "#3A62D9",  
  "AntID_ID_PreThal_SMN" = "#93F6F8", 
  "AntID_ID_PreThal_PMN" = "#808AA9", 
  "AntID_ID_PreThal_MMN" = "#F0A881", 
  "AntID_ID_PreThal_PVH_SON" = "#DFCE93", 
  "MMN_SMN" = "#A700FF", 
  "NPC_Ascl1" = "#00FF1B", 
  "NPC_Neurog2" = "#4EF860"
)

# Filter the color palette to only include the present brain regions
filtered_colors <- full_color_palette[names(full_color_palette) %in% present_brain_regions]

# Plot with activators on the right and repressors on the left
p <- ggraph(ggraph_layout) + 
  geom_edge_link(aes(color = type), 
                 arrow = arrow(length = unit(3, 'mm'), type = "closed"),  # Set a fixed length for all arrows
                 end_cap = circle(3, 'mm'), edge_width = 1) +
  geom_node_point(aes(x = x, y = y, color = BrainRegion), size = 3) +  # Color nodes by brain region
  geom_node_text(aes(x = x, y = y, label = name), repel = TRUE, size = 5, max.overlaps = 200) +  # Label nodes
  scale_edge_color_manual(values = c("Activator" = "blue", "Repressor" = "red")) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_color_manual(values = filtered_colors) +  # Only include colors for present regions
  theme_void() +
  labs(title = "Pathway Plot of Selected eRegulons (Top 10 Percentile by Triplet Rank) with Brain Regions",
       edge_color = "Type",
       edge_width = "Importance TF2G",
       color = "Brain Region")

# Print the plot
print(p)
postscript(file = "Figures/Nkx2-2/SCENIC_Nkx2-2_eRegulon_Pathway_Regulon_TF.eps", width = 30, height = 15,
           horizontal = FALSE, onefile = FALSE, paper = "special")
print(p)
dev.off()

####GRN PLot-old####
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

# Load traG2M NPCription factor list
tf_list_url <- "http://humantfs2.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt"
tf_list <- fread(tf_list_url, header = FALSE, sep = "\t")
traG2M NPCription_factors <- tf_list$V1

# Function to capitalize only the first character
capitalize_first <- function(gene) {
  paste0(toupper(substr(gene, 1, 1)), tolower(substr(gene, 2, nchar(gene))))
}

# Apply the function to the list of genes
formatted_genes <- sapply(traG2M NPCription_factors, capitalize_first)

# Load eRegulon data
data <- read.table("SCENIC/eRegulon/PATTERN_eRegulons_extended.tsv", header = TRUE, sep = "\t")

# Filter and prepare the data
filtered_data <- data %>%
  filter(grepl("extended_\\+/\\+", eRegulon_name) | grepl("extended_\\-/\\+", eRegulon_name)) %>%
  select(eRegulon_name, Gene, importance_TF2G, triplet_rank) %>%
  mutate(type = ifelse(grepl("\\+/\\+", eRegulon_name), "Activator", ifelse(grepl("\\-/\\+", eRegulon_name), "Repressor", NA)))

# Select specific eRegulons
selected_eRegulons <- c("Nkx2-2_extended_+/+", "Nkx2-2_extended_-/+")

filtered_data <- filtered_data %>%
  filter(eRegulon_name %in% selected_eRegulons)

# Combine 'Nkx2-2_extended_+/+' and 'Nkx2-2_extended_-/+' into 'Nkx2-2'
filtered_data <- filtered_data %>%
  mutate(eRegulon_name = str_replace(eRegulon_name, "_extended_\\+/\\+", ""),
         eRegulon_name = str_replace(eRegulon_name, "_extended_\\-/\\+", ""))

# Convert eRegulon_name to a factor with levels ordered by the original eRegulon_name
filtered_data <- filtered_data %>%
  mutate(eRegulon_name = factor(eRegulon_name, levels = unique(eRegulon_name)))

# Keep only traG2M NPCription factors from the formatted list
filtered_data <- filtered_data[filtered_data$Gene %in% formatted_genes, ]

# Normalize the triplet_rank for arrow length mapping
filtered_data <- filtered_data %>%
  mutate(normalized_rank = (max(triplet_rank) - triplet_rank) / max(triplet_rank) * 4 + 1)  # Scale to range 1-5 for arrow length

# Load expression data to extract brain regions
lists <- read.csv(file = "Figures/Original/Fig1/Pattern_scRNA/DEG_pattern_scRNA_compact.csv",
                  header = T, row.names = 1)


desired_clusters <- c("AntID_ID", 
                      "MMN", "Neural Pro (Ascl1)", "Neural Pro (Neurog2)", "NPC", "NPC (Ascl1)", 
                      "NPC (gliogenic)", "NPC (Neurog2)", "G2M NPC",
                      "PMN", "POA/SCN", "PreThal", "PVH_SON", "SMN", "Tub")

# Filter genes with avg_log2FC > 2
filtered_genes <- lists %>%
  filter(cluster %in% desired_clusters) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(avg_log2FC > 2)

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
                                "NPC (gliogenic)" = "#E6B3B3", "NPC (Neurog2)" = "#C3AED3", "G2M NPC" = "#AC8DC4",
                                "G2M NPC (Neurog2)" = "#C3AED3", "G2M NPC_NPC" = "#4EFAA8", 
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
postscript(file = "Figures/Nkx2-2/SCENIC_Nkx2-2_eRegulon_Pathway_Regulon_TF.eps", width = 30, height = 15,
           horizontal = FALSE, onefile = FALSE, paper = "special")
print(p)
dev.off()

####Dot plot####
write.csv(filtered_data, "Figures/Nkx2-2/Nkx2-2_filtered_data.csv", row.names = FALSE)

#rank_threshold <- quantile(filtered_data$triplet_rank, 0.30)

# Filter data to select only the top 10 percentile based on triplet_rank
#filtered_data2 <- filtered_data %>%
#  filter(triplet_rank <= rank_threshold)

#Find top ranked GEnes
ordered_genes <- filtered_data %>%
  select(Gene, type, triplet_rank) %>%  # Include triplet_rank in the selection
  arrange(type, triplet_rank) %>%  # Arrange by type, then by triplet_rank within each type
  distinct(Gene, .keep_all = TRUE) %>%  # Ensure all columns are kept for distinct genes
  pull(Gene)  # Pull the Gene column as a vector

# Print the ordered gene list to verify
print(ordered_genes)

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
meta.data <- Mutant@meta.data
graph <- 100*prop.table(table(Idents(Mutant), Mutant$Genotype), margin = 1) #margin 1= row, 2 = column
dat <- melt(graph)
write.csv(dat, file = "Figures/Nkx2-2/dat.csv", row.names = F)

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

library(scCustomize)
dot_plot <- DotPlot_scCustom(seurat_object = Mutant, 
                             features = ordered_genes, group.by = "Genotype", x_lab_rotate = TRUE, y_lab_rotate = TRUE)
ggsave("Figures/Nkx2-2/SCENIC_Nkx2-2_dotplot.eps", plot = dot_plot, width = 20, height = 8, units = "in")
print(dot_plot)
dev.off()

#2 Heatmap
library(pheatmap)
valid_genes <- ordered_genes[ordered_genes %in% rownames(Mutant[["SCT"]])]
Scaling  <- ScaleData(Mutant, features = valid_genes )
Avg.Scaling <- AverageExpression(Scaling, assay = "SCT", slot = "counts",
                                 verbose = T) #exponential minus 1 mean(expm1(x)) of the SCTTransform value
newnames <- lapply(
  rownames(Avg.Scaling$SCT[valid_genes,]),
  function(x) bquote(italic(.(x))))
heatmp <- pheatmap(Avg.Scaling$SCT[valid_genes,], fontsize = 8, scale = "row",
                   fontsize_col = 8, cluster_rows = F,
                   cluster_cols = F, cellwidth = 10, cellheight = 10, angle_col = 45,
                   color = colorRampPalette(c("blue", "white", "red"))(50),
                   labels_row = as.expression(newnames))

#heatmap v2
expr_data <- Avg.Scaling$SCT[valid_genes, ]
zscore_data <- t(apply(expr_data, 1, function(x) (x - mean(x)) / sd(x)))
heatmp <- pheatmap(zscore_data, fontsize = 8, scale = "none",
                   fontsize_col = 8, cluster_rows = FALSE,
                   cluster_cols = FALSE, cellwidth = 10, cellheight = 10,
                   angle_col = 45, color = colorRampPalette(c("blue", "white", "red"))(50),
                   labels_row = as.expression(newnames))
library(grDevices)
eps_file <- "Figures/Nkx2-2/Nkx2-2_heatmap.eps"
setEPS()
postscript(eps_file, width = 8, height = 15) 
print(heatmp)
dev.off()

#Find unique regions
Mutant <- FindNeighbors(Mutant, dims = 1:20, reduction = "harmony")
Mutant <- FindClusters(Mutant, resolution = 0.8) #SCT_snn_res.0.8
DimPlot(Mutant, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

contingency_table <- table(Mutant@active.ident, Mutant$Genotype)
percentage_table <- prop.table(contingency_table, margin = 1) * 100

library(tidyr)
library(dplyr)
reshaped_df <- as.data.frame(percentage_table) %>%
  pivot_wider(names_from = Var2, values_from = Freq, id_cols = Var1)
reshaped_df$Mutant_Control_Ratio <- with(reshaped_df, Mutant / Control)
View(reshaped_df)

clusters <- c("15","5","9","8","10","6","14")
for (cluster in clusters) {
  # Find markers for the current cluster
  markers <- FindMarkers(Mutant, ident.1 = cluster,  # No need to convert if already a character
                         logfc.threshold = 0.5,
                         min.pct = 0.2, verbose = TRUE)
  
  # Define the file path using %s for string formatting
  file_name <- sprintf("Figures/Nkx2-2/Nkx2-2_Enriched_Mutant_C%s.csv", cluster)
  
  # Save to CSV
  write.csv(markers, file_name, row.names = T)
}

#p <- FeaturePlot_scCustom(seurat_object = Mutant, features = c("Lhx9"),
#                          alpha_exp = 0.75, split.by = "Genotype") & NoLegend() & NoAxes() +
#  theme(plot.title = element_blank())
#png(filename = "Figures/Nkx2-2/Lhx9.png", width = 3000, height = 3000, res = 300)
#print(p)
#dev.off()

plot_and_save_genes <- function(seurat_object, features, directory = "Figures") {
  # Ensure the directory exists
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  for (gene in features) {
    # Generate the plot for the current gene
    p <- FeaturePlot_scCustom(seurat_object = seurat_object, features = c(gene),
                              order = T, split.by = "Genotype", pt.size = 0.1) & NoLegend() & NoAxes() +
      theme(plot.title = element_blank())
    
    # Define the file path
    file_path <- file.path(directory, paste0(gene, ".png"))
    
    # Save the plot to a PNG file
    png(filename = file_path, width = 2400, height = 1200, res = 300)
    print(p)
    dev.off()
  }
}

genes <- c("Dcx","Pmch","Hmx2","Isl1","Dlx1","Dlx2","Prox1","Lef1","Nr4a2","Foxp2","Cartpt",
  "Calb2","Onecut2","Otp","Sim1","Chchd10","Calb1","Fezf1","Foxb1","Lhx1","Fezf2")
Stacked_VlnPlot(seurat_object = Mutant, features = genes, x_lab_rotate = TRUE, group.by = "Cluster")
VlnPlot(Mutant, genes, group.by = "Cluster")

genes <- unique(filtered_data$Gene)
plot_and_save_genes(Mutant, genes, "Figures/Nkx2-2")
