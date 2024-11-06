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
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram
setwd("/media/thomaskim/Data/")

####pattern-matching scATAC####
load(file = "scRNA/Robj/pattern.Robj")

Idents(pattern) <- "Cluster_Pass2"
pattern <- RenameIdents(pattern, "Ant_ID" = "AntID_ID", 
                        "MMN (Prog)" = "MMN", 
                        "PMN (Pro)" = "PMN",
                        "PMN? (Gsx1)" = "PMN",
                        "PreThal (Pro)" = "PreThal", 
                        "PreThal (Sp8)" = "PreThal",
                        "PreThal (Sp9)" = "PreThal",
                        "PreThal_AntID" = "PreThal",
                        "PVH_SON_AntVen" = "PVH_SON",
                        "SMN (Prog)" = "SMN",
                        "Tub (ARC_VMH)" = "Tub", 
                        "Tub (LH-Pmch)"= "Tub",
                        "Tub (LH)"= "Tub",
                        "Tub (Prog)"= "Tub",
                        "NPC (gliogenic?)" = "NPC (gliogenic)",
                        "PMN? (Gsx1)" = "PMN (Gsx1)",
                        "POA/SCN?" = "POA/SCN",
                        "SST (DMH?)" = "SST (DMH)",
                        "ZLI?" = "ZLI")

# Define the mapping between regions and their respective genes
region_gene_map <- list("AntID_ID" = "Lhx6", "MMN" = "Foxb1", "PMN" = "Hmx2",
                        "POA_SCN" = "Cartpt", "PreThal" = "Meis2", "PVH_SON" = "Otp",
                        "SMN" ="Lmx1b", "Tub" = "Nr5a1")

# Function to create directory if it doesn't exist
create_dir_if_not_exists <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# Loop through each region and its corresponding gene
for (region in names(region_gene_map)) {
  gene <- region_gene_map[[region]]
  
  # Generate the plot
  p <- FeaturePlot(pattern, features = gene, order = TRUE, raster = FALSE) + NoLegend() + NoAxes()
  
  # Construct the file name
  filename <- paste0("Figures/Fig1/Pattern_scRNA/Pattern_", region, "_", gene, ".png")
  
  # Create directory if it doesn't exist
  create_dir_if_not_exists(dirname(filename))
  
  # Save the plot as a PNG file
  png(filename = filename, width = 900, height = 900)
  print(p)
  dev.off()
}

# List of genes to plot
genes_to_plot <- c("Sp8","Sp9","Gsx1","Hmx2","Hmx3","Gal", "Cartpt", "Bsx",
                   "Drd2","Ghrh","Th","Gfra2","Aldh1a3","Sst","Npy","Prdm12",
                   "Pnoc","Pmch","Npvf","Hcrt")

# Function to create directory if it doesn't exist
create_dir_if_not_exists <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# Loop through each gene
for (gene in genes_to_plot) {
  
  # Generate the plot
  p <- FeaturePlot(pattern, features = gene, order = TRUE, raster = FALSE) + NoLegend() + NoAxes()
  
  # Construct the file name
  filename <- paste0("Figures/Fig1/Pattern_scRNA/Pattern_", gene, ".png")
  
  # Create directory if it doesn't exist
  create_dir_if_not_exists(dirname(filename))
  
  # Save the plot as a PNG file
  png(filename = filename, width = 900, height = 900)
  print(p)
  dev.off()
}

####neurogenesis-matching scATAC####
load(file = "scRNA/Robj/neurogenesis.Robj")

Idents(neurogenesis) <- "Cluster_Pass2"
neurogenesis <- subset(neurogenesis, idents = c("ARC (Npy_Agrp)", "ARC (Pomc)", "ARC_VMH", "ARC_VMH (Tac2_PMN)", 
                                                "Astrocytes",  "DMH-LH (Grp_Cck)", 
                                                "DMH-LH (Npvf)", "DMH (Npw_PMN)", "Ependymal", "Isl1 cells", 
                                                "LH (Hcrt_Oxt)", "LH (Pmch_Trh)", "MMN", "MMN (Cck)", "MMN (Npy_Cck)", 
                                                "MMN (Nts_Tac1_Cck)", "Neural Pro (Ascl1/Neurog2)", "Neural Pro (Neurog2)", 
                                                "NPC (Epen)", "NPC (Epen/Tan)", "NPC (Glial)", "NPC (Glial_Astro)", 
                                                "NPC (Glial_Check)", "NPC (Glial_Oligo)", "NPC (OPC)", "Oligodendrocytes (newly)", 
                                                "OPC", "PMN (Ghrh_Gal_Cited1)", "PMN (Ghrh_Gal_Oxt)", "POA_SCN", 
                                                "PreThal", "PreThal_ID", "PVN_SON (Avp_Gal_Oxt)", "SCN (Rorb)", 
                                                "SMN", "SMN (Calb2)", "SMN (Tac2)", "Sst_Npy", "Tanycyte", 
                                                "TMN_PH (Hdc)"))
neurogenesis <- RenameIdents(neurogenesis,
                             "NPC (Epen)" = "NPC (Glial)",
                             "NPC (Epen/Tan)" = "NPC (Glial)", 
                             "NPC (Glial)" = "NPC (Glial)",
                             "NPC (Glial_Astro)" = "NPC (Glial)", 
                             "NPC (Glial_Check)" = "NPC (Glial)",
                             "NPC (Glial_Oligo)" = "NPC (Glial)",
                             "NPC (OPC)" = "NPC (Glial)")
neurogenesis <- RenameIdents(neurogenesis,   "NPC (Glial)" = "NPC (Glia)",
                             "ARC (Npy_Agrp)" = "ARC_VMH_PMN",
                             "ARC (Pomc)" = "ARC_VMH_PMN",
                             "ARC_VMH" = "ARC_VMH",
                             "ARC_VMH (Tac2_PMN)" = "ARC_VMH_KISS1",
                             "Astrocytes" = "Astrocyte",
                             "DMH-LH (Grp_Cck)"  = "PMN",
                             "DMH-LH (Npvf)" = "Npvf_Lhx9",
                             "DMH (Npw_PMN)"  = "PMN",
                             "Ependymal" = "Tanycyte_Ependymal",
                             "Isl1 cells" = "ARC_VMH_PMN",
                             "LH (Hcrt_Oxt)" = "Npvf_Lhx9",
                             "LH (Pmch_Trh)" = "ARC_VMH_PMN_PMCH",
                             "MMN" = "MMN",
                             "MMN (Cck)"  = "MMN",
                             "MMN (Npy_Cck)"  = "MMN",
                             "MMN (Nts_Tac1_Cck)"  = "MMN",
                             "Neural Pro (Ascl1/Neurog2)" = "Neural Pro (Neurog2)",
                             "Neural Pro (Neurog2)" = "Neural Pro (Neurog2)",
                             "Oligodendrocytes (newly)" = "Oligodendrocyte",
                             "OPC"  = "OPC",
                             "PMN (Ghrh_Gal_Cited1)" = "PMN",
                             "PMN (Ghrh_Gal_Oxt)" = "PMN",
                             "POA_SCN" = "SCN",
                             "PreThal" = "PreThal",
                             "PreThal_ID" = "PreThal_AntID_ZI",
                             "PVN_SON (Avp_Gal_Oxt)" = "PVH_SON",
                             "SCN (Rorb)" = "SCN",
                             "SMN" = "SMN",
                             "SMN (Calb2)" = "SMN",
                             "SMN (Tac2)" = "SMN",
                             "Sst_Npy" = "ARC_VMH_PMN",
                             "Tanycyte" = "Tanycyte_Ependymal",
                             "TMN_PH (Hdc)" = "ARC_VMH_PMN" )

# Define the mapping between regions and their respective genes
region_gene_map <- list("NPC (Glia)" = "Slc7a10",
                        "ARC1" = "Kiss1",
                        "ARC2" = "Agrp",
                        "ARC3" = "Pmch",
                        "MMN" = "Nts",
                        "LH" = "Npvf",
                        "PMN1" = "Ghrh",
                        "PreThal" = "Meis2",
                        "PVH" = "Avp",
                        "SCN" = "Penk",
                        "SMN" = "Tac1")

# Function to create directory if it doesn't exist
create_dir_if_not_exists <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# Loop through each region and its corresponding gene
for (region in names(region_gene_map)) {
  gene <- region_gene_map[[region]]
  
  # Generate the plot
  p <- FeaturePlot(neurogenesis, features = gene, order = F, raster = FALSE) + NoLegend() + NoAxes()
  
  # Construct the file name
  filename <- paste0("Figures/Fig1/neurogenesis_scRNA/neurogenesis_", region, "_", gene, ".png")
  
  # Create directory if it doesn't exist
  create_dir_if_not_exists(dirname(filename))
  
  # Save the plot as a PNG file
  png(filename = filename, width = 900, height = 900)
  print(p)
  dev.off()
}

# List of genes to plot
genes_to_plot <- c("Agrp","Avp","Kiss1","Tac2","Npvf","Th","Prph","Npw","Pmch",
                   "Sst","Npy","Pomc","Ghrh","Hdc","Nts","Cck","Gal","Pnoc","Penk",
                   "Hcrt","Tac1","Tac2","Grp","Oxt","Trh","Crh","Vip","Cartpt","Ar","Ghr","Prlr","Bsx")

# Function to create directory if it doesn't exist
create_dir_if_not_exists <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# Loop through each gene
for (gene in genes_to_plot) {
  
  # Generate the plot
  p <- FeaturePlot(neurogenesis, features = gene, order = TRUE, raster = FALSE) + NoLegend() + NoAxes()
  
  # Construct the file name
  filename <- paste0("Figures/Fig1/neurogenesis_scRNA/neurogenesis_", gene, ".png")
  
  # Create directory if it doesn't exist
  create_dir_if_not_exists(dirname(filename))
  
  # Save the plot as a PNG file
  png(filename = filename, width = 900, height = 900)
  print(p)
  dev.off()
}
