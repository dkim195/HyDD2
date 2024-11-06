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

####pattern####
load(file = "scRNA/Robj/pattern.Robj")

####pattern-matching scATAC####
Idents(pattern) <- "Cluster_Pass2"
pattern <- subset(pattern, idents = c("NPC", "NPC (Ascl1)", "NPC (Neurog2)"))
pattern <- SCTransform(pattern,
                       vars.to.regress = c("nCount_RNA","nFeature_RNA"))
pattern <- RunPCA(pattern, npcs = 50, ndims.print = NA, verbose = F)
pattern <- RunHarmony(pattern, group.by.vars = "orig.ident",
                      assay.use = "SCT", plot.convergence = T)
pattern <- RunUMAP(pattern, reduction = "harmony", dims = 1:10,
                   n.neighbours = 30L, min.dist = 0.1, spread = 1)
DimPlot(pattern, reduction = "umap", label = T,pt.size = 0.1) + NoLegend() + NoAxes()

Idents(pattern) <- "Age_Sum"
pattern <- RenameIdents(pattern, "E11" = "Early", "E12" = "Early", "E13" = "Late", "E14" = "Late")
pattern <- AddMetaData(pattern, pattern@active.ident, "NPCAge")
pattern$CompositeLabel <- paste(pattern$Cluster_Pass2, pattern$NPCAge, sep = "_")
Idents(pattern) <- "CompositeLabel"
identities <- Idents(pattern)
new_order <- c("NPC_Early","NPC_Late","NPC (Neurog2)_Early",
               "NPC (Neurog2)_Late","NPC (Ascl1)_Early","NPC (Ascl1)_Late")
Idents(pattern) <- factor(identities, levels = new_order)

#find DEG
markers <- FindAllMarkers(pattern, test.use = "wilcox",
                          logfc.threshold = 1,
                          min.pct = 0.1, verbose = T)
write.csv(markers, file = "Figures/DEG_NPC.csv")


library(pheatmap)
#pattern <- RenameIdents(pattern, "E11" = "Early", "E12" = "Early", "E13" = "Late", "E14" = "Late")
genes <- c("Lhx6", "Prdm12", "Onecut3", "Foxd1", 
           "Gad1",  "Sp9", "Dlx6", 
           "Lhx1", "Lhx2","Arx", "Lhx1os", "Slc32a1", 
           "Foxb1",  "Nkx2-2", "Nkx2-1","Uncx", "Sim1", 
           "Rax", "Ascl1", "Cdkn2c", "Ccnb1", "Ube2c", 
           "Gadd45g", "Wnt9a", "Wnt8b", "Neurod4", "Neurog1", 
           "Neurog2", "Gsx1", "Btg2", "Ppp1r14a", "Nhlh1", "Olig1", "Neurod1", 
           "Neurod2", 
           "Hmx2", "Hmx3",  "Cited1", "Bsx", "Calb1",  "Prox1", "Cartpt",  "Gbx2", 
           "Otx1", "Zic1", "Lef1", "Zic4", "Otp", "Insm2", "Arxes2", 
           "Meis2", "Sp8",  "Gad2", 
           "Isl1", "Dlx1", "Dlx6", "Dlx2","Dlx5", "Esrrg", "Irx6", "Nr4a2", "Lmx1b", 
           "Irx5", "Irx3", "Barhl1", "Lmx1a", "Foxa1", 
           "Neurod6", "Nr5a1", "Sox14", "Nr0b1", "Chchd10", 
           "Prdm13", "Six6", "Lin28a", 
           "Hmga1", "Hmgn2", "Crabp2", "Igdcc3", "Nfia", "Nfib", "Nfix", "Sox8", "Sox9", "Hes1", "Id1", "Id4",
           "Notch1", "Zbtb20", "Fabp7", "Rfx4", "Tcf4", "Tox3", "Npas3",
           "Sox6", "Tsc22d4", "Nr2f1", "Nr2f2",
           "Rspo3", "Zic4", "Zic5", "Fgf15", "Pax6","Fezf1", "Pcp4", "Nkx2-3", "Foxa2","Sp5","Foxa1","Barhl1","Irx3","Lmx1b","Lmx1a")
genes <- unique(genes)
Scaling  <- ScaleData(pattern, features = genes)
Avg.Scaling <- AverageExpression(Scaling, assay = "SCT", slot = "counts",
                                 verbose = T) #exponential minus 1 mean(expm1(x)) of the SCTTransform value
newnames <- lapply(
  rownames(Avg.Scaling$SCT[genes,]),
  function(x) bquote(italic(.(x))))

heatmp <- pheatmap(Avg.Scaling$SCT[genes,], fontsize = 8, scale = "row",
                   fontsize_col = 8, cluster_rows = T,
                   cluster_cols = T, cellwidth = 10, cellheight = 7, angle_col = 45,
                   color = colorRampPalette(c("blue", "white", "red"))(50),
                   labels_row = as.expression(newnames))
heatmp
library(grDevices)
eps_file <- "Figures/Fig1/heatmap_NPC.eps"
setEPS()
postscript(eps_file, width = 8, height = 15) 
print(heatmp)
dev.off()



####SCENIC####