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

# Set random seed for reproducibility
set.seed(1234)

# Configure parallelization
plan("multicore", workers = 10)
future.seed = TRUE

# Increase maximum allowed memory size
options(future.globals.maxSize = 512000 * 1024^2) # 512 GB RAM

# Set working directory
setwd("/faststorage/project/Hypothalamus/data")

# Load the pre-saved Seurat object
load(file = "scATAC/Isl1_Mutant_adjusted_Final.Robj")
DefaultAssay(Mutant_adjusted) <- 'peaks'
pfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE",
                                            tax_group = 'vertebrates', all_versions = FALSE))
Mutant_adjusted <- AddMotifs(object = Mutant_adjusted,
                             genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm, assay = 'peaks')
Mutant_adjusted <- RunChromVAR(
  object = Mutant_adjusted,
  genome = BSgenome.Mmusculus.UCSC.mm10)

# Save the updated Seurat object
save(Mutant_adjusted, file = "scATAC/Isl1_Mutant_adjusted_Final2.Robj")

# Set default assay to 'chromvar'
DefaultAssay(Mutant_adjusted) <- 'chromvar'

# Define region to motif mapping
region_motif_map <- list(
  "REGION1" = "MA0668.1", "REGION2" = "MA0890.1", "REGION3" = "MA0669.1", 
  "REGION4" = "MA0720.1", "REGION5" = "MA0701.2", "REGION6" = "MA0802.1",
  "REGION7" = "MA0678.1", "REGION8" = "MA0631.1", "REGION9" = "MA0800.1", "REGION10" = "MA0898.1", "REGION11" = "MA0897.1" )

            

# Function to create directory if it doesn't exist
create_dir_if_not_exists <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# Create the main directory if it doesn't exist
output_dir <- "/faststorage/project/Hypothalamus/plots/Fig1/"
create_dir_if_not_exists(output_dir)

# Loop through each region and its corresponding motif
for (region in names(region_motif_map)) {
  motif <- region_motif_map[[region]]
  
  # Generate the plot
  p <- FeaturePlot(Mutant_adjusted, features = motif, order = TRUE, raster = FALSE, 
                   min.cutoff = 'q10', max.cutoff = 'q90', split.by = "Genotype")
  
  # Construct the file name
  filename <- paste0(output_dir, "Isl1_", region, "_", motif, ".png")
  
  # Save the plot as a PNG file
  png(filename = filename, width = 2400, height = 1200)
  print(p)
  dev.off()
}

# Define the list of motifs for the MotifPlot
MOTIFNAMES <- c("MA0668.1", "MA0890.1", "MA0669.1", "MA0720.1", "MA0701.2", "MA0802.1", "MA0678.1", "MA0631.1", "MA0800.1", "MA0898.1", "MA0897.1")

# Generate the motif plot
p <- MotifPlot(
  object = Mutant_adjusted,
  motifs = MOTIFNAMES,
  assay = 'peaks'
)

# Save the motif plot as a PNG file
filename <- paste0(output_dir, "Isl1_MOTIFLIST.png")
png(filename = filename, width = 3000, height = 3000, res = 300)
print(p)
dev.off()
