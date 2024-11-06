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
library(data.table)
library(igraph)
library(ggraph)
library(ggrepel)
library(pheatmap)
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram
setwd("/media/thomaskim/Data/")

#Modified from TEST1 and TEST2

####Function####
# Function to find feedback pairs between two motifs
FoundFeedBackPairsOne <- function(tmp_motif1, tmp_motif2){
  ##### Positive and Negative pairs from tmp_motif1 ###############
  tmp_motif1_pos <- tmp_motif1[tmp_motif1$type == 'Activator', ]
  tmp_motif1_neg <- tmp_motif1[tmp_motif1$type == 'Repressor', ]
  
  ##### Positive and Negative pairs from tmp_motif2 ###############
  tmp_motif2_pos <- tmp_motif2[tmp_motif2$type == 'Activator', ]
  tmp_motif2_neg <- tmp_motif2[tmp_motif2$type == 'Repressor', ]
  
  ##### Find overlapping positive pairs ##########
  Pos_index1_rev <- paste(tmp_motif1_pos$Gene, tmp_motif1_pos$eRegulon_name)
  Pos_index1_rev_index <- which(Pos_index1_rev %in% paste(tmp_motif2_pos$eRegulon_name, tmp_motif2_pos$Gene))
  
  if(length(Pos_index1_rev_index) > 0){
    Pos_index1_rev_res <- tmp_motif1_pos[Pos_index1_rev_index, ]
    Pos_index1_rev_res <- Pos_index1_rev_res[!duplicated(Pos_index1_rev_res$Gene), ]
    Pos_index1_rev_res <- Pos_index1_rev_res[, c("eRegulon_name", "Gene", "type")]
  } else {
    Pos_index1_rev_res <- data.frame(eRegulon_name = 'ND', Gene = 'ND', type = 'ND')
  }
  
  ##### Find overlapping negative pairs ##########
  Neg_index1_rev <- paste(tmp_motif1_neg$Gene, tmp_motif1_neg$eRegulon_name)
  Neg_index1_rev_index <- which(Neg_index1_rev %in% paste(tmp_motif2_neg$eRegulon_name, tmp_motif2_neg$Gene))
  
  if(length(Neg_index1_rev_index) > 0){
    Neg_index1_rev_res <- tmp_motif1_neg[Neg_index1_rev_index, ]
    Neg_index1_rev_res <- Neg_index1_rev_res[!duplicated(Neg_index1_rev_res$Gene), ]
    Neg_index1_rev_res <- Neg_index1_rev_res[, c("eRegulon_name", "Gene", "type")]
  } else {
    Neg_index1_rev_res <- data.frame(eRegulon_name = 'ND', Gene = 'ND', type = 'ND')
  }
  
  #### Combine Positive and Negative pairs ##########
  Res <- rbind(Pos_index1_rev_res, Neg_index1_rev_res)
  Res_nd <- which(Res$eRegulon_name == 'ND')
  if(length(Res_nd) > 0){
    Res <- Res[-Res_nd, ]
  }
  
  return(Res)
}

# Function to find feedback pairs across all motifs
FoundFeedBackPairs_new <- function(Motif_list){
  Out_list <- list()
  for(i in 1:length(Motif_list)){
    for(j in 1:length(Motif_list)){
      print(c(i, j))
      Names_1 <- names(Motif_list)[i]
      Names_2 <- names(Motif_list)[j]
      print(c(Names_1, Names_2))
      tmp_motif1 <- Motif_list[[i]]
      tmp_motif2 <- Motif_list[[j]]
      tmp_out <- FoundFeedBackPairsOne(tmp_motif1, tmp_motif2)
      
      # Check if tmp_out has rows before adding columns
      if (nrow(tmp_out) > 0) {
        tmp_out$celltypes <- paste(Names_1, Names_2, sep = ':')
        tmp_out$TFs <- paste(Names_1, tmp_out$eRegulon_name, sep = ':')
        tmp_out$Target <- paste(Names_2, tmp_out$Gene, sep = ':')
        Out_list <- c(Out_list, list(tmp_out))
      }
    }
  }
  Out_list_out <- do.call("rbind", Out_list)
  print(dim(Out_list_out))
  return(Out_list_out)
}

# Annotate feedback results
Process_the_Feedback_res <- function(RPCMG_Feedback_res){
  sp_TFs <- strsplit(RPCMG_Feedback_res$TFs, split = ':', fixed = TRUE)
  sp_Target <- strsplit(RPCMG_Feedback_res$Target, split = ':', fixed = TRUE)
  RPCMG_Feedback_res$TF_index <- sapply(sp_TFs, function(x) x[[1]])
  RPCMG_Feedback_res$TF_gene <- sapply(sp_TFs, function(x) x[[2]])
  RPCMG_Feedback_res$Target_index <- sapply(sp_Target, function(x) x[[1]])
  RPCMG_Feedback_res$Target_gene <- sapply(sp_Target, function(x) x[[2]])
  RPCMG_Feedback_res$TF_target <- paste(RPCMG_Feedback_res$TF_gene, RPCMG_Feedback_res$Target_gene, sep = '::')
  return(RPCMG_Feedback_res)
}

####Select genes####
# Keep only unique genes in each group
NPC_Ascl1_Early <- c("Meis2","Arx","Sp8","Ascl1","Sp9","Dlx5","Hmx3","Sox14","Pax6","Nr5a1","Nkx2-2","Nr2f1")
NPC_Ascl1_Late <- c("Neurod1","Lhx6","Prox1","Lhx1","Otp","Nkx2-1","Neurod6")
NPC_Neurog2_Early <- c("Irx5","Irx3","Foxa2","Sp5","Nhlh1", "Sim1", "Otx1")
NPC_Neurog2_Late <- c("Sox9","Rfx4","Nfia","Nfib","Sox6","Nfix","Tcf4")
NPC_Early <- c("Nr0b1","Zic4","Zic5","Zic3","Zic1","Nr2e1","Fezf1")
NPC_Late <-c("Lhx2","Rax","Pcp4","Sox8","Onecut3","Six6")


genes <- c("Lhx6", "Prdm12", "Onecut3", "Foxd1", 
           "Gad1",  "Sp9", "Dlx6", 
           "Lhx1", "Lhx2","Arx", "Lhx1os", "Slc32a1", 
           "Foxb1",  "Nkx2-2", "Nkx2-1","Uncx", "Sim1", 
           "Rax", "Ascl1", "Cdkn2c", "Ccnb1", "Ube2c", 
           "Gadd45g", "Wnt9a", "Wnt8b", "Neurod4", "Neurog1", 
           "Neurog2", "Gsx1", "Btg2", "Ppp1r14a", "Nhlh1", "Olig1", "Neurod1", 
           "Neurod2", 
           "Hmx2", "Hmx3",  "Cited1", "Bsx", "Calb1",  "Prox1", "Cartpt",  "Gbx2", 
           "Otx1", "Zic1", "Lef1", "Zic4", "Otp", "Insm2", "Arxes2", 
           "Meis2", "Sp8",  "Gad2", 
           "Isl1", "Dlx1", "Dlx6", "Dlx2","Dlx5", "Esrrg", "Irx6", "Nr4a2", "Lmx1b", 
           "Irx5", "Irx3", "Barhl1", "Lmx1a", "Foxa1", 
           "Neurod6", "Nr5a1", "Sox14", "Nr0b1", "Chchd10", 
           "Prdm13", "Six6", "Lin28a", 
           "Hmga1", "Hmgn2", "Crabp2", "Igdcc3", "Nfia", "Nfib", "Nfix", "Sox8", "Sox9", "Hes1", "Id1", "Id4",
           "Notch1", "Zbtb20", "Fabp7", "Rfx4", "Tcf4", "Tox3", "Npas3",
           "Sox6", "Tsc22d4", "Nr2f1", "Nr2f2",
           "Rspo3", "Zic4", "Zic5", "Fgf15", "Pax6","Fezf1", "Pcp4", "Nkx2-3", "Foxa2","Sp5","Foxa1","Barhl1","Irx3","Lmx1b","Lmx1a")


# Combine all unique genes
All_genes_test_unique <- c(NPC_Ascl1_Early, NPC_Ascl1_Late, NPC_Neurog2_Early,
                           NPC_Neurog2_Late, NPC_Early, NPC_Late)

#filter eREgulon
data <- read.table("SCENIC/eRegulon/PATTERN_eRegulons_extended.tsv", header = TRUE, sep = "\t")
filtered_data <- data %>%
  filter(grepl("extended_\\+/\\+", eRegulon_name) | grepl("extended_\\-/\\+", eRegulon_name)) %>%
  select(eRegulon_name, Gene, importance_TF2G, triplet_rank) %>%
  mutate(type = ifelse(grepl("\\+/\\+", eRegulon_name), "Activator", ifelse(grepl("\\-/\\+", eRegulon_name), "Repressor", NA)))
filtered_data <- filtered_data %>%
  mutate(eRegulon_name = str_replace(eRegulon_name, "_extended_\\+/\\+", ""),
         eRegulon_name = str_replace(eRegulon_name, "_extended_\\-/\\+", ""))

f1 <- which(filtered_data$eRegulon_name %in% NPC_Ascl1_Early == TRUE)
f2 <- which(filtered_data$eRegulon_name %in% NPC_Ascl1_Late == TRUE)
f3 <- which(filtered_data$eRegulon_name %in% NPC_Neurog2_Early == TRUE)
f4 <- which(filtered_data$eRegulon_name %in% NPC_Neurog2_Late == TRUE)
f5 <- which(filtered_data$eRegulon_name %in% NPC_Early == TRUE)
f6 <- which(filtered_data$eRegulon_name %in% NPC_Late  == TRUE)


NPC_Ascl1_Early_Reg_motif <- filtered_data[f1, ]
NPC_Ascl1_Late_Reg_motif <- filtered_data[f2, ]
NPC_Neurog2_Early_Reg_motif <- filtered_data[f3, ]
NPC_Neurog2_Late_Reg_motif <- filtered_data[f4, ]
NPC_Early_Reg_motif <- filtered_data[f5, ]
NPC_Late_Reg_motif <- filtered_data[f6, ]

#GRN
GRNs_list <- list(NPC_Ascl1_Early_Reg_motif, NPC_Ascl1_Late_Reg_motif, NPC_Neurog2_Early_Reg_motif,
                  NPC_Neurog2_Late_Reg_motif, NPC_Early_Reg_motif, NPC_Late_Reg_motif)
names(GRNs_list) <- c("NPC_(Ascl1)_Early", "NPC_(Ascl1)_Late", "NPC_(Neurog2)_Early",
                      "NPC_(Neurog2)_Late", "NPC_Early", "NPC_Late")
GRNs_list_cl <- GRNs_list

#forming GRN
Early_Feedback_res <- FoundFeedBackPairs_new(GRNs_list_cl)
# Sort the feedback pairs by type if needed (optional step)
Early_Feedback_res <- Early_Feedback_res[order(Early_Feedback_res$type), ]
# Annotate the feedback pairs
Early_Feedback_res <- Process_the_Feedback_res(Early_Feedback_res)

#PLOT
graph_data <- Early_Feedback_res %>%
  select(TFs, Target, type) %>%
  distinct()
# Create an igraph object
graph <- graph_from_data_frame(graph_data, directed = TRUE)
# Plot the graph using ggraph
ggraph(graph, layout = 'fr') +  # You can choose different layouts like 'fr', 'circle', etc.
  geom_edge_link(aes(color = type), arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() +
  scale_edge_color_manual(values = c("Activator" = "blue", "Repressor" = "red")) +
  labs(title = "Feedback TF Pairs Network")


#plot 2
graph_data <- Early_Feedback_res %>%
  mutate(source_region = sapply(strsplit(TFs, ":", fixed = TRUE), function(x) x[1]),
         target_region = sapply(strsplit(Target, ":", fixed = TRUE), function(x) x[1])) %>%
  select(source_region, target_region, type) %>%
  distinct()
# Create an igraph object from the simplified data
graph <- graph_from_data_frame(graph_data, directed = TRUE)
# Plot the graph using ggraph
ggraph(graph, layout = 'fr') +  # You can choose different layouts like 'fr', 'circle', etc.
  geom_edge_link(aes(color = type), arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() +
  scale_edge_color_manual(values = c("Activator" = "blue", "Repressor" = "red")) +
  labs(title = "Region Interaction Network: Activation and Repression")

####Test####
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

# Step 1: Add Cluster labels to both dataframes
NPC_Ascl1_Early_Reg_motif$Cluster <- "NPC_(Ascl1)_Early"
NPC_Ascl1_Late_Reg_motif$Cluster <- "NPC_(Ascl1)_Late"
NPC_Neurog2_Early_Reg_motif$Cluster <- "NPC_(Neurog2)_Early"
NPC_Neurog2_Late_Reg_motif$Cluster <- "NPC_(Neurog2)_Late"
NPC_Early_Reg_motif$Cluster <- "NPC_Early"
NPC_Late_Reg_motif$Cluster <- "NPC_Late"

# Step 2: Combine the dataframes into one
combined_motif <- bind_rows(NPC_Ascl1_Early_Reg_motif, NPC_Ascl1_Late_Reg_motif, NPC_Neurog2_Early_Reg_motif,
                            NPC_Neurog2_Late_Reg_motif, NPC_Early_Reg_motif, NPC_Late_Reg_motif)

# Keep only transcription factors from the formatted list
combined_motif <- combined_motif[combined_motif$Gene %in% formatted_genes, ]

# Step 3: Identify conflicting genes within each cluster separately
conflicting_genes <- combined_motif %>%
  group_by(Cluster, Gene) %>%
  filter(n_distinct(type) > 1) %>%
  pull(Gene) %>%
  unique()

# Step 4: Exclude conflicting genes within each cluster
cleaned_combined_motif <- combined_motif %>%
  filter(!(Gene %in% conflicting_genes & Cluster == Cluster))  # Ensure filtering within the same cluster

# Step 5: Modify the source column to represent the Cluster
cleaned_combined_motif$source <- cleaned_combined_motif$Cluster

# Step 6: Create and plot the combined regulatory network
# Create an igraph object from the cleaned data
graph_combined_cleaned <- graph_from_data_frame(cleaned_combined_motif %>% select(source, Gene, type, importance_TF2G), directed = TRUE)

# Plotting the graph using ggraph
ggraph(graph_combined_cleaned, layout = 'fr') + 
  geom_edge_link(aes(color = type), 
                 arrow = arrow(length = unit(cleaned_combined_motif$importance_TF2G / max(cleaned_combined_motif$importance_TF2G) * 4 + 1, 'mm'), type = "closed"), 
                 end_cap = circle(3, 'mm'))  +  
  geom_node_text(aes(label = name), repel = TRUE, max.overlaps = 100) +  # Increased max.overlaps
  scale_edge_color_manual(values = c("Activator" = "blue", "Repressor" = "red")) +
  theme_void() +
  labs(title = "Combined Regulatory Network for SMN and MMN Clusters",
       edge_color = "Type",
       edge_width = "Importance TF2G")