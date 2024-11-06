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

#load(file = "scATAC/neurogenesis_atac_annotation_v3b.Robj")

#neurogenesis_atac <- RunChromVAR(
#  object = neurogenesis_atac,
#  genome = BSgenome.Mmusculus.UCSC.mm10)

#save(neurogenesis_atac, file = "scATAC/neurogenesis_atac_annotation_v4.Robj")

load(file = "scATAC/neurogenesis_atac_annotation_v4.Robj")

DefaultAssay(neurogenesis_atac) <- 'chromvar'
# Define the mapping between regions and their respective genes
#region_motif_map <- list("NPC (Glia)" = "MA0671.1",
#                         "ARC1" = "MA1529.1",
#                         "ARC2" = "MA0592.3",
#                         "ARC3" = "MA0711.1",
#                         "MMN" = "MA0756.1",
#                         "LH" = "MA1577.1",
#                         "PMN1" = "MA1608.1",
#                         "PreThal" = "MA1476.1",
#                         "PVH" = "MA0887.1",
#                         "SCN" = "MA1151.1",
#                         "SMN" = "MA0702.2")

region_motif_map <- list("NPC (Glia)" = "MA0161.2",
                         "ARC1" = "MA0643.1",
                         "ARC2" = "MA1608.1",
                         "ARC3" = "MA0648.1",
                         "MMN" = "MA0756.1",
                         "LH" = "MA0075.3",
                         "PMN1" = "MA0893.2",
                         "PreThal" = "MA0881.1",
                         "PVH" = "MA0793.1",
                         "SCN" = "MA1150.1",
                         "SMN" = "MA0703.2")


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
  p <- FeaturePlot(neurogenesis_atac, features = motif, order = TRUE, raster = FALSE, 
                   min.cutoff = 'q10', max.cutoff = 'q90') + NoLegend() + NoAxes()
  
  # Construct the file name
  filename <- paste0("/faststorage/project/Hypothalamus/plots/Fig1/neurogenesis_scATAC/neurogenesis_", region, "_", motif, ".png")
  
  # Create directory if it doesn't exist
  create_dir_if_not_exists(dirname(filename))
  
  # Save the plot as a PNG file
  png(filename = filename, width = 900, height = 900)
  print(p)
  dev.off()
}

MOTIFNAMES <- c("MA0671.1","MA1529.1","MA0592.3","MA0711.1","MA0756.1","MA1577.1","MA1608.1","MA1476.1","MA0887.1","MA1151.1","MA0702.2",
                "MA0161.2","MA0643.1","MA1608.1","MA0648.1","MA0756.1","MA0075.3","MA0893.2","MA0881.1","MA0793.1","MA1150.1","MA0703.2")

p <- MotifPlot(
  object = neurogenesis_atac,
  motifs = MOTIFNAMES,
  assay = 'peaks'
)

filename <- paste0("/faststorage/project/Hypothalamus/plots/Fig1/neurogenesis_scATAC/neurogenesis_MOTIFLIST.png")
png(filename = filename, width = 3000, height = 3000, res = 300)
print(p)
dev.off()