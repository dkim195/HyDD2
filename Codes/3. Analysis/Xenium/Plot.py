import os
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
os.chdir("/media/thomaskim/Data/Xenium")

#P21_HET_REGION2
adata = sc.read_h5ad(filename="OUTPUT/P21_HOMO_REGION4.h5ad")

#Plot4
#Plot4
genes_to_plot = ["Atoh7", "Olig2", "Olig1",  "Gjc3", "Opalin", "Glul", "Ntsr2",
  "Gfap", "Hcrt", "Pvalb", "Gad2", "Meis2", 
  "Kiss1", "Calb2", "Ar", "Cartpt", "Hdc", 
  "Lhx1", "Lhx6", "Nkx2-2", "Nr5a1", "Sim1", 
  "Slc17a6", "Nr2f2", "Pax6", "Pnoc", "Prdm13", 
  "Penk", "Trh", "Fos"]


for gene in genes_to_plot:
    # Create a figure with specified size
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plotting function with adjusted parameters for each gene
    sq.pl.spatial_scatter(
        adata, library_id="spatial",
        color=[gene], shape=None,
        size=0.1, cmap="Reds", img=True, img_alpha=0.6, ax=ax,
        legend_loc="right", legend_fontsize="xx-small"
    )

    # Rotate x-axis labels for better readability
    for label in ax.get_xticklabels():
        label.set_rotation(90)

    # Adjust layout
    plt.tight_layout()

    # Save the figure with the gene name in the filename
    filename = f"/media/thomaskim/Data/Figures/Xenium/P21_HOMO_REGION4_{gene}.png"
    fig.savefig(filename, format='png', dpi=300, bbox_inches='tight')

    # Display the plot
    plt.show()

    # Close the figure to free up memory
    plt.close(fig)

