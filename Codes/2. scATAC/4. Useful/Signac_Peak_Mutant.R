library(cowplot)
library(dplyr)
library(tibble)
library(Matrix)
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
set.seed(1234)
setwd("/media/thomaskim/Data/scATAC/")

#### Directory Setup ####
ensure_dir_exists <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}
# Base directory for saving outputs
base_dir <- "/media/thomaskim/Data/scATAC/CSV/"
experiment_name <- "Lhx2_E12_scATAC"  # Adjust as necessary
# Full path construction
full_path <- file.path(base_dir, experiment_name)
# Ensure the directory exists
ensure_dir_exists(full_path)


#### Load Data ####
load(file = "Files/Lhx2_Mutant_adjusted_Final.Robj")

#### Motif Analysis ####
# Load necessary libraries for motif analysis
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
# Retrieve motifs and annotate
pfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE",
                                            tax_group = 'vertebrates', all_versions = FALSE))
Mutant_adjusted <- AddMotifs(object = Mutant_adjusted,
                   genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm, assay = 'peaks')

DefaultAssay(Mutant_adjusted) <- 'peakunion'
my_clusters <- unique(Mutant_adjusted$Cluster_Pass2)
Idents(Mutant_adjusted) <- "Cluster_Pass2"

# Map function to perform both DE analysis and adding gene names
cluster_markers <- purrr::map(my_clusters, function(cluster) {
  # Perform differential expression analysis
  markers <- FindMarkers(Mutant_adjusted, 
                         ident.1 = "Control", 
                         ident.2 = "Mutant",
                         group.by = "Genotype", 
                         subset.ident = cluster,
                         logfc.threshold = 0.1, min.pct = 0.05, 
                         test.use = "wilcox")
  
  # Temporarily change DefaultAssay to 'peaks' for ClosestFeature
  DefaultAssay(Mutant_adjusted) <- "peaks"
  
  # Proceed only if markers are found
  if (nrow(markers) > 0) {
    # Get gene names for the rownames in the marker data frame
    Genename <- ClosestFeature(Mutant_adjusted, regions = rownames(markers))
    # Add gene names to markers data frame
    markers$gene_name <- Genename$gene_name
  }
  
  # Restore DefaultAssay back to 'peakunion' for the next iteration or subsequent analysis
  DefaultAssay(Mutant_adjusted) <- 'peakunion'
  
  return(markers)
})
# Optionally, reset DefaultAssay to a default state if needed outside the loop
DefaultAssay(Mutant_adjusted) <- 'peakunion'
# Naming each list element with cluster IDs if not already done
names(cluster_markers) <- my_clusters

# Loop through each DE result
for (i in seq_along(cluster_markers)) {
  # Process each DE result with dplyr for cleaner syntax
  cluster_markers[[i]] <- cluster_markers[[i]] %>%
    dplyr::relocate(gene_name, .before = 1) %>%
    {
      if("p_val" %in% names(.)) dplyr::select(., -p_val) else .
    } %>%
    dplyr::rename(., Control = pct.1, Mutant = pct.2) %>%
    dplyr::arrange(., avg_log2FC)
  
  # Construct file path with actual cluster name
  file_name <- sprintf("DE_result_%s.csv", names(cluster_markers)[i])
  file_path <- file.path(full_path, file_name)
  
  # Save to CSV
  write.csv(cluster_markers[[i]], file_path, row.names = FALSE)
}

####MOTIF####
#Mutant_adjusted <- RegionStats(Mutant_adjusted, genome = BSgenome.Mmusculus.UCSC.mm10)
# Processing each cluster
DefaultAssay(Mutant_adjusted) <- 'peaks'
for (i in seq_along(cluster_markers)) {
  # Extract current cluster markers
  current_markers <- cluster_markers[[i]]
  
  # Filter for control-enriched peaks (positive avg_log2FC)
  control_enriched_peaks <- rownames(current_markers[current_markers$avg_log2FC > 2, ])
  
  # Filter for mutant-enriched peaks (negative avg_log2FC)
  mutant_enriched_peaks <- rownames(current_markers[current_markers$avg_log2FC < -2, ])
  
  # Conduct motif enrichment analysis for control-enriched peaks
  if (length(control_enriched_peaks) > 0) {
    enriched_motifs_control <- FindMotifs(object = Mutant_adjusted, features = control_enriched_peaks)
    # Save or process enriched_motifs_control
    motifs_control_file_path <- file.path(full_path, sprintf("Control_Enriched_Motifs_Cluster_%s.csv", names(cluster_markers)[i]))
    write.csv(as.data.frame(enriched_motifs_control), motifs_control_file_path, row.names = FALSE)
  }
  
  # Conduct motif enrichment analysis for mutant-enriched peaks
  if (length(mutant_enriched_peaks) > 0) {
    enriched_motifs_mutant <- FindMotifs(object = Mutant_adjusted, features = mutant_enriched_peaks)
    # Save or process enriched_motifs_mutant
    motifs_mutant_file_path <- file.path(full_path, sprintf("Mutant_Enriched_Motifs_Cluster_%s.csv", names(cluster_markers)[i]))
    write.csv(as.data.frame(enriched_motifs_mutant), motifs_mutant_file_path, row.names = FALSE)
  }
}





MotifPlot(object = Mutant_adjusted, motifs = "MA1489.1")
table(Mutant_adjusted$Genotype, Mutant_adjusted$Cluster_Pass2)
FeaturePlot(Mutant_adjusted, "Foxb1", split.by = "Genotype",   max.cutoff = 'q95')

