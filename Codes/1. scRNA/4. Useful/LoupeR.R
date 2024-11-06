library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(loupeR)
set.seed(1234)
setwd("/media/thomaskim/Data/")

# Gene Expression RNA assay
load(file = "scRNA/Robj/Nkx21_Final.Robj")
Mutant <- UpdateSeuratObject(Mutant)
Mutant@meta.data$old.ident <- NULL
Mutant@meta.data$orig.ident <- NULL
Mutant@meta.data$Cluster_Pass1 <- NULL
Mutant@meta.data$Reference <- NULL
Mutant@meta.data$predicted.id <- NULL
Idents(Mutant) <- "Genotype"

#get assay
count_mat <- Mutant@assays$RNA@counts
clusters <- select_clusters(Mutant)
projections <- select_projections(Mutant)


# convert the count matrix, clusters, and projections into a Loupe file
create_loupe(count_mat = count_mat, 
             clusters = clusters, 
             projections = projections, 
             output_name = "scRNA/create_loupe_test")
