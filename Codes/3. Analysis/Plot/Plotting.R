####Plotting Split Genes####
library(ggplot2)
p <- FeaturePlot(Mutant, features = c("Foxb1"), split.by = "Genotype",
                 pt.size = 0.5, order = F, combine = F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()  + ggtitle("")  + theme(legend.position = "none")
}
cowplot::plot_grid(plotlist = p)

####Stacked Violin####
library(ggplot2)
####test
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0),
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot.margin )
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45, hjust = 0.5), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(y)) +
                            expand_limits(y = y))
  
  patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
}

StackedVlnPlot(Mutant, c("Foxb1","Lhx1","Lhx9") ) +
  theme(axis.text.x = element_text(angle = 45, hjust=))

####Volcano####

lists <- read.csv(file = "Figures/Original/Fig1/Pattern_scRNA/DEG_pattern_scRNA_compact.csv",
                  header = T, row.names = 1)
desired_clusters <- c("AntID_ID", 
                      "MMN")

# Filter genes with avg_log2FC > 2
filtered_genes <- lists %>%
  filter(cluster %in% desired_clusters) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) #%>%
# filter(avg_log2FC > 1)

ggplot(data = filtered_genes, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                  label=gene))  + geom_point() + ggtitle("NRCAM") +
  geom_text(hjust=-1,vjust=1)


genes.to.label1 = c("Pax6","Dlx1")

ggplot(data = filtered_genes, aes(x = avg_log2FC, y = -log10(p_val_adj),label=gene))  + geom_point() + ggtitle("NRCAM") +
  geom_text_repel(data=subset(filtered_genes, gene %in% c(genes.to.label1)),
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label=gene),
                  color="red", size = 10, nudge_y      = 0.05,
                  vjust        = 0,
                  segment.size = 1.2) +
  geom_point(data=subset(filtered_genes, gene %in% c(genes.to.label1)),
             aes(x = avg_log2FC, y = -log10(p_val_adj)), color="red",size=3) +
  labs(x = "avg_log2FC", y = "-log10(p_val_adj)")

