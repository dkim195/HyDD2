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
library(scCustomize)
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram
setwd("/media/thomaskim/Data/")


#single gene
count_gene_expression_by_genotype <- function(seurat_obj, gene_of_interest, genotype_column = "Genotype",
                                              control_label = "Control", mutant_label = "Mutant", expression_threshold = 0) {
  
  # 1. Check if the gene exists in the dataset
  if (!(gene_of_interest %in% rownames(seurat_obj))) {
    stop(paste("Gene", gene_of_interest, "not found in dataset."))
  }
  
  # 2. Check if the genotype column exists in the meta.data
  if (!(genotype_column %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Column", genotype_column, "not found in metadata."))
  }
  
  # 3. Fetch the expression data for the gene of interest
  expression_data <- FetchData(seurat_obj, vars = gene_of_interest, slot = "counts")
  
  # 4. Subset the data into Control and Mutant groups
  control_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[genotype_column]] == control_label, ])
  mutant_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[genotype_column]] == mutant_label, ])
  
  # 5. Identify cells expressing the gene above the threshold in each group
  control_expr_cells <- control_cells[expression_data[control_cells, gene_of_interest] > expression_threshold]
  mutant_expr_cells <- mutant_cells[expression_data[mutant_cells, gene_of_interest] > expression_threshold]
  
  # 6. Count the number of cells expressing the gene in each Genotype
  control_expr_count <- length(control_expr_cells)
  mutant_expr_count <- length(mutant_expr_cells)
  
  # 7. Return the counts as a named list
  return(list(
    Control = control_expr_count,
    Mutant = mutant_expr_count
  ))
}

# Example usage:
# result <- count_gene_expression_by_genotype(Mutant, "Isl1", genotype_column = "Genotype",control_label = "Control", mutant_label = "Mutant", expression_threshold = 1)
# print(result)



#multiple genes
count_gene_expression_by_genotype_multiple_genes <- function(seurat_obj, genes_of_interest,
                                                             genotype_column = "Genotype",
                                                             control_label = "Control", mutant_label = "Mutant",
                                                             expression_threshold = 0, output_file = "gene_expression_counts.csv") {
  
  # Initialize an empty list to store the results
  result_list <- list()
  
  for (gene_of_interest in genes_of_interest) {
    
    # 1. Check if the gene exists in the dataset
    if (!(gene_of_interest %in% rownames(seurat_obj))) {
      warning(paste("Gene", gene_of_interest, "not found in dataset. Skipping..."))
      next
    }
    
    # 2. Fetch the expression data for the gene of interest
    expression_data <- FetchData(seurat_obj, vars = gene_of_interest, slot = "counts")
    
    # 3. Subset the data into Control and Mutant groups
    control_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[genotype_column]] == control_label, ])
    mutant_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[genotype_column]] == mutant_label, ])
    
    # 4. Identify cells expressing the gene above the threshold in each group
    control_expr_cells <- control_cells[expression_data[control_cells, gene_of_interest] > expression_threshold]
    mutant_expr_cells <- mutant_cells[expression_data[mutant_cells, gene_of_interest] > expression_threshold]
    
    # 5. Count the number of cells expressing the gene in each Genotype
    control_expr_count <- length(control_expr_cells)
    mutant_expr_count <- length(mutant_expr_cells)
    
    # 6. Add the counts to the result list
    result_list[[gene_of_interest]] <- data.frame(
      Gene = gene_of_interest,
      Control = control_expr_count,
      Mutant = mutant_expr_count
    )
  }
  
  # Combine the individual data frames into a single data frame
  result_df <- do.call(rbind, result_list)
  
  # 7. Save the results to a CSV file
  write.csv(result_df, file = output_file, row.names = FALSE)
  
  # Return the data frame as the function output
  return(result_df)
}

# Example usage:
# genes_to_check <- c("Isl1", "Th", "Npy")
# result <- count_gene_expression_by_genotype_multiple_genes(Mutant, genes_to_check, genotype_column = "Genotype", control_label = "Control", mutant_label = "Mutant", expression_threshold = 1, output_file = "gene_expression_results.csv")
# print(result)

genes_to_check <- c("Isl1", "Th", "Npy")
result <- count_gene_expression_by_genotype_multiple_genes(Mutant, genes_to_check,
                                                           genotype_column = "Genotype", control_label = "Control",
                                                           mutant_label = "Mutant", expression_threshold = 1,
                                                           output_file = "gene_expression_results.csv")
print(result)


####Isl1####
load(file = "scRNA/Robj/Isl1_Final.Robj")
VlnPlot(Mutant, "Isl1", group.by = "Genotype", add.noise = F, sort = T,
        same.y.lims = T, log = T, pt.size = 0) + NoLegend() + NoAxes()

genes_to_check <- c("Isl1", "Tbr1", "Slc17a6")
result <- count_gene_expression_by_genotype_multiple_genes(Mutant, genes_to_check,
                                                           genotype_column = "Genotype", control_label = "Control",
                                                           mutant_label = "Mutant", expression_threshold = 0,
                                                           output_file = "Figures/cKO value/Isl1_gene_expression_results.csv")
print(result)


####Lhx1####
load(file = "scRNA/Robj/Lhx1_Final.Robj")
VlnPlot(Mutant, "Lhx1", group.by = "Genotype", add.noise = F, sort = T,
        same.y.lims = T, log = T, pt.size = 0) + NoLegend() + NoAxes()

genes_to_check <- c("Isl1", "Ebf2", "Gad1", "Neurog2", "Ascl1","Foxb1",
                    "Sim1", "Foxa2", "Lmx1b", "Barhl2", "Calb2", "Shox2", "Nr4a2", "Tbr1")
result <- count_gene_expression_by_genotype_multiple_genes(Mutant, genes_to_check,
                                                           genotype_column = "Genotype", control_label = "Control",
                                                           mutant_label = "Mutant", expression_threshold = 0,
                                                           output_file = "Figures/cKO value/Lhx1_gene_expression_results.csv")
print(result)

plot_and_save_genes <- function(seurat_object, features, directory = "Figures") {
  # Ensure the directory exists
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  for (gene in features) {
    # Generate the plot for the current gene
    p <- FeaturePlot_scCustom(seurat_object = seurat_object, features = c(gene),
                              order = T, split.by = "Genotype", pt.size = 0.1) & NoLegend() & NoAxes() +
      theme(plot.title = element_blank())
    
    # Define the file path
    file_path <- file.path(directory, paste0(gene, ".png"))
    
    # Save the plot to a PNG file
    png(filename = file_path, width = 2400, height = 1200, res = 300)
    print(p)
    dev.off()
  }
}

genes <-  c("Neurog2","Ascl1")
plot_and_save_genes(Mutant, genes, "Figures/Lhx1")

####Nkx2-2####
load(file = "scRNA/Robj/Nkx22_Final.Robj")
VlnPlot(Mutant, "Nkx2-2", group.by = "Genotype", add.noise = F, sort = T,
        same.y.lims = T, log = T, pt.size = 0) + NoLegend() + NoAxes()

genes_to_check <- c("Nkx2-2","Foxb1","Hmx2","Otp", "Zic1", "Irx5", "Sp9", "Pitx2", "Sp8","Meis2","Dlx1",
                    "Arx", "Pitx2", "Barhl1", "Lmx1a", "Pmch", "Sim1", "Otp", "Lhx1", "Dlx2")
result <- count_gene_expression_by_genotype_multiple_genes(Mutant, genes_to_check,
                                                           genotype_column = "Genotype", control_label = "Control",
                                                           mutant_label = "Mutant", expression_threshold = 0,
                                                           output_file = "Figures/cKO value/Nkx22_gene_expression_results.csv")


plot_and_save_genes <- function(seurat_object, features, directory = "Figures") {
  # Ensure the directory exists
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  for (gene in features) {
    # Generate the plot for the current gene
    p <- FeaturePlot_scCustom(seurat_object = seurat_object, features = c(gene),
                              order = T, split.by = "Genotype", pt.size = 0.1) & NoLegend() & NoAxes() +
      theme(plot.title = element_blank())
    
    # Define the file path
    file_path <- file.path(directory, paste0(gene, ".png"))
    
    # Save the plot to a PNG file
    png(filename = file_path, width = 2400, height = 1200, res = 300)
    print(p)
    dev.off()
  }
}

genes <-  c("Foxb1","Hmx2","Otp", "Zic1", "Irx5", "Sp9", "Pitx2", "Sp8","Meis2")
plot_and_save_genes(Mutant, genes, "Figures/Nkx2-2/")


####Dlx1/2####
load(file = "scRNA/Robj/Dlx1_2_Final.Robj")
VlnPlot(Mutant, "Dlx1", group.by = "Genotype", add.noise = F, sort = T,
        same.y.lims = T, log = T, pt.size = 0) + NoLegend() + NoAxes()
VlnPlot(Mutant, "Dlx2", group.by = "Genotype", add.noise = F, sort = T,
        same.y.lims = T, log = T, pt.size = 0) + NoLegend() + NoAxes()

genes_to_check <- c("Dlx1","Dlx2","Slc32a1","Nr4a2","Arx","Gad2","Isl1","Dlx5",
                    "Foxb1","Irx5","Sox14","Foxa1","Hmx2","Meis2","Sp8","Sp9","Pitx2","Pax6")
result <- count_gene_expression_by_genotype_multiple_genes(Mutant, genes_to_check,
                                                           genotype_column = "Genotype", control_label = "Control",
                                                           mutant_label = "Mutant", expression_threshold = 0,
                                                           output_file = "Figures/cKO value/Dlx_gene_expression_results.csv")

plot_and_save_genes <- function(seurat_object, features, directory = "Figures") {
  # Ensure the directory exists
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  for (gene in features) {
    # Generate the plot for the current gene
    p <- FeaturePlot_scCustom(seurat_object = seurat_object, features = c(gene),
                              order = T, split.by = "Genotype", pt.size = 0.1) & NoLegend() & NoAxes() +
      theme(plot.title = element_blank())
    
    # Define the file path
    file_path <- file.path(directory, paste0(gene, ".png"))
    
    # Save the plot to a PNG file
    png(filename = file_path, width = 2400, height = 1200, res = 300)
    print(p)
    dev.off()
  }
}

genes <-   c("Dlx1","Dlx2","Slc32a1","Nr4a2","Arx","Gad2","Isl1","Dlx5",
             "Foxb1","Irx5","Sox14","Foxa1","Hmx2","Meis2","Sp8","Sp9","Pitx2","Pax6")
plot_and_save_genes(Mutant, genes, "Figures/Dlx/")



####P8 Dlx1/2####
load(file = "scRNA/Robj/P8_Dlx_Final.Robj")
Idents(Mutant) <- "Cluster_Pass3"
#modifying based on 1. location and neuropeptide
#Annotation from Mutant P8 is not super great
Mutant <- RenameIdents(Mutant, "ARC (Npw_Agrp)" = "ARC (Npy_Agrp)", 
                       "ARC (Npy_Agrp)" = "ARC (Npy_Agrp)", 
                       "ARC (Pomc)" = "ARC (Pomc)", 
                       "DMH (Npw_PMN)" = "DMH_PMN (Npw)",
                       "DMH-LH (Grp_Cck)" = "DMH_PMN (Grp_Bsx)", 
                       "LH (Hcrt)" = "LH_DMH (Hcrt_Npvf)",
                       "LH (Pmch_Tac2)" = "ARC_LH (Kiss1_Tac2_Pmch)", 
                       "LH (Pmch_Trh)" = "PVH_SON (Oxt_Trh_Crh)",
                       "MMN (Calb1_Cck)" = "MMN", 
                       "MMN (Cartpt_Cck)" = "MMN" ,
                       "MMN (Cck)" = "MMN", 
                       "MMN (Nts_Cartpt)" = "MMN",
                       "MMN (Nts_Tac1_Cck)" = "MMN", 
                       "MMN (Pcp4_Cck)" = "MMN",
                       "PMN (Ghrh_Gal_Oxt)" = "ARC_PMN (Ghrh_Gal_Th)" ,
                       "PreThal" = "PreThal",
                       "PVN_SON (Avp_Gal_Oxt)" = "PVN_SON (Avp)", 
                       "SCN (Avp)" = "PVN_SON (Avp)", 
                       "SCN (Rorb)" = "SCN (Nms)",
                       "SCN (Vip)" = "SCN (Vip)",
                       "SMN (Calb2_Tac1)" = "SMN",
                       "SMN (Calb2)" = "Thalamus", 
                       "SMN (Nts)" = "SMN",
                       "Sst (ARC-TMN?)" = "UNKNOWN (Sst_Pthlh)" ,
                       "Sst (TMN-PH?)" = "TMN_PH (Sst_Thrb)",
                       "Sst_Npy"  = "ARC (Sst_Npy_Th)",
                       "Tac2" = "SMN",
                       "TMN_PH (Hdc)" = "TMN_PH (Hdc_Prph)",
                       "VMH"  = "VMH (Nr5a1)",
                       "VMH (Tac1)"  = "VMH (Nr5a1_Tac1)")
identities <- Idents(Mutant)
new_order <- sort(levels(identities))
Idents(Mutant) <- factor(identities, levels = new_order)

VlnPlot(Mutant, "Dlx1", group.by = "Genotype", add.noise = F, sort = T,
        same.y.lims = T, log = T, pt.size = 0) + NoLegend() + NoAxes()
VlnPlot(Mutant, "Dlx2", group.by = "Genotype", add.noise = F, sort = T,
        same.y.lims = T, log = T, pt.size = 0) + NoLegend() + NoAxes()

check <- table(Mutant@active.ident, Mutant$Genotype)
write.csv(check, file = "Figures/cKO value/P8_Dlx.csv")

genes_to_check <- c("Dlx1","Dlx2","Slc32a1","Nr4a2","Arx","Gad2","Isl1","Dlx5",
                    "Foxb1","Irx5","Sox14","Foxa1","Hmx2","Meis2","Sp8","Sp9","Pitx2","Pax6",
                    "Vip","Nms","Fgf8","Prok2","Vipr2","Six6","Rora","Rorb","Slc17a6","Cck","Tac1","Slc17a7","Pnoc","Nts",
                    "Npy","Cartpt","Npw","Bsx","Sst","Agrp","Tac2","Shox2","Gbx2","Gad1","Th","Gal","Ghrh","Qrfp","Lhx9","Hcrt","Npvf")
result <- count_gene_expression_by_genotype_multiple_genes(Mutant, genes_to_check,
                                                           genotype_column = "Genotype", control_label = "Control",
                                                           mutant_label = "Mutant", expression_threshold = 0,
                                                           output_file = "Figures/cKO value/P8Dlx_gene_expression_results.csv")

plot_and_save_genes <- function(seurat_object, features, directory = "Figures") {
  # Ensure the directory exists
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  for (gene in features) {
    # Generate the plot for the current gene
    p <- FeaturePlot_scCustom(seurat_object = seurat_object, features = c(gene),
                              order = T, split.by = "Genotype", pt.size = 0.1) & NoLegend() & NoAxes() +
      theme(plot.title = element_blank())
    
    # Define the file path
    file_path <- file.path(directory, paste0(gene, ".png"))
    
    # Save the plot to a PNG file
    png(filename = file_path, width = 1200, height = 600, res = 100)
    print(p)
    dev.off()
  }
}

genes <- c("Dlx1","Dlx2","Slc32a1","Nr4a2","Arx","Gad2","Isl1","Dlx5",
             "Foxb1","Irx5","Sox14","Foxa1","Hmx2","Meis2","Sp8","Sp9","Pitx2","Pax6",
             "Vip","Nms","Fgf8","Prok2","Vipr2","Six6","Rora","Rorb","Slc17a6","Cck","Tac1","Slc17a7","Pnoc","Nts",
             "Npy","Cartpt","Npw","Bsx","Sst","Agrp","Tac2","Shox2","Gbx2","Gad1","Th","Gal","Ghrh","Qrfp","Lhx9","Hcrt","Npvf",
           "Pthlh")
plot_and_save_genes(Mutant, genes, "Figures/P8Dlx/")
