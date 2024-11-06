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
library(cetcolor)
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram
setwd("/media/thomaskim/Data/")

####Dlx_scRNA####
#plot
load(file = "scRNA/Robj/P8_Dlx_Final.Robj")

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

clusters <-c("1","26","12","10","15","4","2","22","5","28","27")
for (cluster in clusters) {
  # Find markers for the current cluster
  markers <- FindMarkers(Mutant, ident.1 = cluster,  # No need to convert if already a character
                         logfc.threshold = 0.5,
                         min.pct = 0.2, verbose = TRUE)
  
  # Define the file path using %s for string formatting
  file_name <- sprintf("Figures/P8Dlx/P8Dlx_Enriched_Mutant_C%s.csv", cluster)
  
  # Save to CSV
  write.csv(markers, file_name, row.names = T)
}

#p <- FeaturePlot_scCustom(seurat_object = Mutant, features = c("Lhx9"),
#                          alpha_exp = 0.75, split.by = "Genotype") & NoLegend() & NoAxes() +
#  theme(plot.title = element_blank())
#png(filename = "Figures/Dlx/Lhx9.png", width = 3000, height = 3000, res = 300)
#print(p)
#dev.off()

library(scCustomize)
plot_and_save_genes <- function(seurat_object, features, directory = "Figures") {
  # Ensure the directory exists
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  for (gene in features) {
    # Generate the plot for the current gene
    p <- FeaturePlot_scCustom(seurat_object = seurat_object, features = c(gene),
                              order = F, split.by = "Genotype", pt.size = 0.1) & NoLegend() & NoAxes() +
      theme(plot.title = element_blank())
    
    # Define the file path
    file_path <- file.path(directory, paste0(gene, ".png"))
    
    # Save the plot to a PNG file
    png(filename = file_path, width = 2400, height = 1200, res = 300)
    print(p)
    dev.off()
  }
}

genes <- c("Npw","Npvf","Hcrt","Lhx9","Qrfp","Ghrh","Gal","Th","Gad1","Gad2","Gbx2","Shox2","Tac2","Agrp","Npy","Sst","Pthlh","Bsx")
Stacked_VlnPlot(seurat_object = Mutant, features = genes, x_lab_rotate = TRUE, group.by = "Cluster")

VlnPlot(Mutant, genes, group.by = "Cluster")
plot_and_save_genes(Mutant, genes, "Figures/P8Dlx")
