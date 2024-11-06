library(cowplot)
library(dplyr)
library(Matrix)
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(future)
library(harmony)
library(presto)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)
plan("multicore", workers = 10)
future.seed = TRUE  # Ensure reproducible random number generation across sessions
options(future.globals.maxSize = 512000 * 1024^2) #512 gb ram
setwd("/faststorage/project/Hypothalamus/data")

#load(file = "scATAC/pattern_atac_annotation_v3b.Robj")

#pattern_atac <- RunChromVAR(
#  object = pattern_atac,
#  genome = BSgenome.Mmusculus.UCSC.mm10)

#save(pattern_atac, file = "scATAC/pattern_atac_annotation_v4.Robj")

load(file = "scATAC/pattern_atac_annotation_v4.Robj")

DefaultAssay(pattern_atac) <- 'chromvar'
# Define the mapping between regions and their respective genes
#region_motif_map <- list("AntID_ID" = "MA0756.1", "MMN" = "MA1529.1", "PMN" = "MA1518.1",
#                         "PreThal" = "MA1476.1", "PVN_SON" = "MA0751.1",
#                        "SMN" ="MA0702.2", "Tub" = "MA0643.1")

#region_motif_map <- list("AntID_ID" = "MA0756.1", "MMN" = "MA0877.2", "PMN" = "MA0889.1",
#                         "PreThal" = "MA0774.1", "PVN_SON" = "MA0696.1",
#                        "SMN" ="MA1112.2", "Tub" = "MA0830.2", "Tub2" = "MA1540.1")


region_motif_map <- list("Lhx6" = "MA0658.1", "Lmx1b" = "MA0703.2", "Foxb1" = "MA0845.1",
                         "Nr5a1" = "MA1540.1", "Meis2" = "MA0774.1",
                        "Hmx2" ="MA0897.1", "Isl1" = "MA1608.1")


# Function to create directory if it doesn't exist
create_dir_if_not_exists <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# Loop through each region and its corresponding motif
for (region in names(region_motif_map)) {
  motif <- region_motif_map[[region]]
  
  # Generate the plot
  p <- FeaturePlot(pattern_atac, features = motif, order = TRUE, raster = FALSE, 
                   min.cutoff = 'q10', max.cutoff = 'q90') + NoLegend() + NoAxes()
  
  # Construct the file name
  filename <- paste0("/faststorage/project/Hypothalamus/plots/Fig1/Pattern_scATAC/Pattern_", region, "_", motif, ".png")
  
  # Create directory if it doesn't exist
  create_dir_if_not_exists(dirname(filename))
  
  # Save the plot as a PNG file
  png(filename = filename, width = 900, height = 900)
  print(p)
  dev.off()
}

#MOTIFNAMES <- c("MA0756.1", "MA1529.1", "MA1518.1",
#                "MA1476.1", "MA0751.1", "MA0702.2", "MA0643.1", "MA0756.1", "MA0877.2", "MA0889.1",
#                "MA0774.1", "MA0696.1", "MA1112.2", "MA0830.2", "MA1540.1")

MOTIFNAMES <- c("MA0658.1", "MA0703.2", "MA0845.1","MA1540.1","MA0774.1", "MA0897.1", "MA1608.1")

p <- MotifPlot(
  object = pattern_atac,
  motifs = MOTIFNAMES,
  assay = 'peaks'
)

filename <- paste0("/faststorage/project/Hypothalamus/plots/Fig1/Pattern_scATAC/Pattern_MOTIFLIST.png")
png(filename = filename, width = 3000, height = 3000, res = 300)
print(p)
dev.off()