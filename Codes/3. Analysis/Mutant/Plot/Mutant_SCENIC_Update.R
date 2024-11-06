####GRN PLot####
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

####TF list####
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

####Load####
# Load eRegulon data
data <- read.table("SCENIC/eRegulon/PATTERN_eRegulons_extended.tsv", header = TRUE, sep = "\t")

# Filter and prepare the data
filtered_data <- data %>%
  filter(grepl("extended_\\+/\\+", eRegulon_name) | grepl("extended_\\-/\\+", eRegulon_name)) %>%
  select(eRegulon_name, Gene, importance_TF2G, triplet_rank) %>%
  mutate(type = ifelse(grepl("\\+/\\+", eRegulon_name), "Activator", ifelse(grepl("\\-/\\+", eRegulon_name), "Repressor", NA)))

####Change####
# Select specific eRegulons
selected_eRegulons <- c("Isl1_extended_+/+", "Isl1_extended_-/+")

####Proceed####
filtered_data <- filtered_data %>%
  filter(eRegulon_name %in% selected_eRegulons)

# Combine 'Isl1_extended_+/+' and 'Isl1_extended_-/+' into 'Isl1'
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

####Define####
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

####Manual Definition based on above#####
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

####Continue####
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
  if (V(graph_eRegulon)$type[i] == "Activator") {
    ggraph_layout$x[i] <- abs(ggraph_layout$x[i])  # Place activators on the right side
  } else if (V(graph_eRegulon)$type[i] == "Repressor") {
    ggraph_layout$x[i] <- -abs(ggraph_layout$x[i])  # Place repressors on the left side
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

