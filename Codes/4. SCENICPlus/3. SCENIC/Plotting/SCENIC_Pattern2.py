#load
import os
import pycisTopic
import pandas as pd
import scanpy as sc
import anndata
import mudata
import matplotlib.pyplot as plt
os.chdir("/faststorage/project/Hypothalamus/SCENIC/PATTERN_SCENIC_outs")
scplus_mdata = mudata.read("PATTERN_scplusmdata.h5mu")

eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,)
eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]

eRegulon_gene_AUC.obs['scRNA_counts:Cluster_Pass2'].cat.categories
unique_groups = eRegulon_gene_AUC.obs['scRNA_counts:Cluster_Pass2'].unique()
scplus_mdata.uns["direct_e_regulon_metadata"]
scplus_mdata.uns["extended_e_regulon_metadata"]

genes = ['Dlx1_direct_+/+_(150g)',
         'Dlx2_direct_+/+_(89g)',
         'Isl1_direct_+/+_(302g)',
         'Lhx1_direct_+/+_(313g)',
         'Lhx2_direct_+/+_(113g)',
         'Nkx2-1_direct_+/+_(117g)',
         'Nkx2-2_direct_+/+_(41g)']

#plotting
# Create a new figure with specific size
fig, ax = plt.subplots(figsize=(10, 6))  # Adjust the size as needed
# Create a plot with specified size directly in Scanpy to handle multiple axes
axes_dict = sc.pl.stacked_violin(eRegulon_gene_AUC, var_names=genes,
                                 groupby='scRNA_counts:Cluster_Pass2', dendrogram=False,
                                 figsize=(10, 6), show=False)

# Rotate the x-axis labels for each axis
for ax in axes_dict.values():
    for label in ax.get_xticklabels():
        label.set_rotation(90)

# Since we have multiple axes, we need to handle the figure from any of the axes
fig = next(iter(axes_dict.values())).get_figure()

# Adjust layout and save
plt.tight_layout()
fig.savefig("pattern_scenic_violin.eps", format='eps', bbox_inches='tight')

# Show the plot
plt.show()

# Close the plot to release resources
plt.close(fig)



#plotting flipped
# Create a new figure with specific size
fig, ax = plt.subplots(figsize=(10, 6))  # Adjust the size as needed
# Create and directly display the plot with axes swapped
axes_dict = sc.pl.stacked_violin(eRegulon_gene_AUC, var_names=genes,
                                 groupby='scRNA_counts:Cluster_Pass2', dendrogram=False,
                                 figsize=(10, 6), swap_axes=True, show=False)

# Since swap_axes=True might change the handling of axes, check if it's still a dictionary
if isinstance(axes_dict, dict):
    # Rotate the x-axis labels for each axis if there are multiple axes
    for ax in axes_dict.values():
        for label in ax.get_xticklabels():
            label.set_rotation(90)
    # Use any axis to get the figure object
    fig = next(iter(axes_dict.values())).get_figure()
else:
    # If only one axis, handle it directly
    for label in axes_dict.get_xticklabels():
        label.set_rotation(90)
    fig = axes_dict.get_figure()

# Adjust layout and save the figure in EPS format
plt.tight_layout()
fig.savefig("pattern_scenic_violin2.eps", format='eps', bbox_inches='tight')

# Show the plot
plt.show()

# Close the plot to release resources
plt.close(fig)















