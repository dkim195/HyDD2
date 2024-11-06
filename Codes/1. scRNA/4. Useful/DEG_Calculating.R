library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(purrr)
set.seed(1234)
setwd("/media/thomaskim/Data/")

####Load####
load(file = "scRNA/Robj/Ctnnb1_Final.Robj")
Mutant <- UpdateSeuratObject(Mutant)

####Directory####
#place to save
ensure_dir_exists <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

base_dir <- "/media/thomaskim/Data/scRNA/CSV/"
seurat_object_name <- "Ctnnb1_E12"

# Full path construction
full_path <- file.path(base_dir, seurat_object_name)

# Ensure the directory exists
ensure_dir_exists(full_path)

####DEG####
# Assuming 'id' is your cluster identifier and 'condition' is your control/mutant label
DimPlot(Mutant, group.by = "Cluster_Pass2", label = T) + NoLegend()

my_clusters <- unique(Mutant$Cluster_Pass2)
Idents(Mutant) <- "Cluster_Pass2"

# Use 'purrr::map' to iterate over each cluster to find markers
cluster_markers <- purrr::map(my_clusters, function(cluster) {
  FindMarkers(Mutant, 
              ident.1 = "Control", 
              ident.2 = "Mutant",
              group.by = "Genotype", 
              subset.ident = cluster,
              logfc.threshold = 0.1, min.pct = 0.05,
              test.use = "wilcox")
})

####Save####
names(cluster_markers) <- my_clusters

# Loop through each DE result
for (i in seq_along(cluster_markers)) {
  # Process each DE result with dplyr for cleaner syntax
  cluster_markers[[i]] <- cluster_markers[[i]] %>%
    # Convert row names (gene names) to a column if not already present
    tibble::rownames_to_column(var = "Gene") %>%
    # Conditionally remove 'p_val' if it exists
    {
      if("p_val" %in% names(.)) select(., -p_val) else .
    } %>%
    # Rename 'pct.1' and 'pct.2'
    rename(Control = pct.1, Mutant = pct.2) %>%
    # Sort by 'avg_log2FC' in ascending order
    arrange(avg_log2FC)
  
  # Construct file path with actual cluster name
  file_name <- sprintf("DE_result_%s.csv", names(cluster_markers)[i])
  file_path <- file.path(full_path, file_name)
  
  # Save to CSV
  write.csv(cluster_markers[[i]], file_path, row.names = FALSE)
}


VlnPlot(Mutant, "Irx2",
        group.by = "Cluster_Pass2", split.by = "Genotype")

####all clusters####
Markers <- FindMarkers(Mutant, 
                       ident.1 = "Control", 
                       ident.2 = "Mutant",
                       group.by = "Genotype", 
                       logfc.threshold = 0.05, min.pct = 0.05,
                       test.use = "wilcox")

marker_name <- "Ctnnb1_Markers" # Replace this with a variable or fixed name as required

# Process the data
processed_data <- Markers %>%
  # Convert row names (gene names) to a column if not already present
  tibble::rownames_to_column(var = "Gene") %>%
  # Conditionally remove 'p_val' if it exists
  {
    if("p_val" %in% names(.)) {
      select(., -p_val)
    } else {
      .
    }
  } %>%
  # Rename 'pct.1' and 'pct.2'
  rename(Control = pct.1, Mutant = pct.2) %>%
  # Sort by 'avg_log2FC' in ascending order
  arrange(avg_log2FC)

# Ensure the directory exists
if (!dir.exists(full_path)) {
  dir.create(full_path, recursive = TRUE)
}

# Create the full file path by appending the marker name and ".csv" to the base path
file_path <- file.path(full_path, paste0(marker_name, ".csv"))

# Save the processed data to this path
write.csv(processed_data, file_path, row.names = FALSE)


####Mutant specific cluster####
Mutant <- FindNeighbors(Mutant, dims = 1:30, reduction = "harmony")
Mutant <- FindClusters(Mutant, resolution = 0.4) #SCT_snn_res.0.8
DimPlot(Mutant, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
table(Mutant$Genotype, Mutant@active.ident)
markers <- FindMarkers(Mutant, ident.1 = "",test.use = "wilcox",
                       logfc.threshold = 0.5,
                       min.pct = 0.1, verbose = T)