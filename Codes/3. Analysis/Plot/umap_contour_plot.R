library(cetcolor)
library(Seurat)
library(ggplot2)

# generate UMAP plot
pl1 <- DimPlot(Mutant, group.by = "Cluster_Pass2", combine = F, split.by = "Genotype")
# custom color scale
scale.col <- cet_pal(16, name = "fire")
# make plot
pl1[[1]] & 
  stat_density_2d(aes_string(x = "UMAP_1", y = "UMAP_2", fill = "after_stat(level)"), 
                  linewidth = 0.2, geom = "density_2d_filled", 
                  colour = "ivory", alpha = 0.4, n = 150, h = c(1.2, 1.2)) & 
  scale_fill_gradientn(colours = scale.col) & 
  DarkTheme()

####thresholding####
# Generate UMAP plot
pl1 <- DimPlot(Mutant, combine = F, cols = group, split.by = "Genotype")
# Custom color scale
scale.col <- cet_pal(16, name = "fire")
# Function to threshold density values
threshold_density <- function(x, from = 0, to = 100) {
  x <- scales::rescale(x, to = c(from, to))
  x[x < from] <- from
  x[x > to] <- to
  return(x)
}
# Make plot with thresholded density values
pl1[[1]] & 
  stat_density_2d(aes_string(x = "UMAP_1", y = "UMAP_2", fill = "threshold_density(after_stat(level))"), 
                  linewidth = 0.2, geom = "density_2d_filled", 
                  colour = "ivory", alpha = 0.4, n = 150, h = c(1.2, 1.2)) & 
  scale_fill_gradientn(colours = scale.col, limits = c(0, 100)) & 
  DarkTheme()
