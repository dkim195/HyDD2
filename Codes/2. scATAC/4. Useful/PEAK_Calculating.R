library(cowplot)
library(dplyr)
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
experiment_name <- "TK98A_Experiment"  # Adjust as necessary

# Full path construction
full_path <- file.path(base_dir, experiment_name)

# Ensure the directory exists
ensure_dir_exists(full_path)

#### Load Data ####
load(file = "Files/TK98A.Robj")

#### Motif Analysis ####
# Load necessary libraries for motif analysis
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)

# Retrieve motifs and annotate
pfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE",
                                            tax_group = 'vertebrates', all_versions = FALSE))

DefaultAssay(TK98A) <- 'peakunion'

TK98A <- AddMotifs(object = TK98A,
                   genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)

#Turn this into Loop
da_peaks <- FindAllMarkers(object = TK98A, 
                           test.use = 'wilcox', logfc.threshold = 0.1, min.pct = 0.05, only.pos = T)
Genename <- ClosestFeature(TK98A, regions = rownames(da_peaks))
da_peaks$gene_name <- Genename$gene_name
# Save differential accessibility results
da_peaks_file_path <- file.path(full_path, "DifferentialAccessibility.csv")
write.csv(da_peaks, da_peaks_file_path, row.names = FALSE)

# peak each cluster MOTIF
#TK98A <- RegionStats(TK98A, genome = BSgenome.Mmusculus.UCSC.mm10)
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$avg_log2FC > 2 & da_peaks$pct.1 > 0.2 & da_peaks$cluster == "8",])

enriched_motifs <- FindMotifs(object = TK98A,
                              features = rownames(da_peaks))
enriched_motifs_file_path <- file.path(full_path, "enriched_motifs.csv")
write.csv(enriched_motifs, enriched_motifs_file_path, row.names = FALSE)


MotifPlot(object = TK98A, motifs = head(rownames(enriched_motifs)))


