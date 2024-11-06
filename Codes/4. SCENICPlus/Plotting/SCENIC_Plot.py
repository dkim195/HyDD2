#load
import os
import pycisTopic
import pandas as pd
import scanpy as sc
import anndata
import mudata
import matplotlib.pyplot as plt
os.chdir("Downloads/SCENIC")
scplus_mdata = mudata.read("scplusmdata.h5mu")

#check
print(scplus_mdata)
print(scplus_mdata['scRNA_counts'].obs)
print(scplus_mdata['scRNA_counts'].var)  
print(scplus_mdata.uns)
eRegulon_gene_AUC.obs['scRNA_counts:Clusters'].cat.categories
unique_groups = eRegulon_gene_AUC.obs['scRNA_counts:Clusters'].unique()
scplus_mdata.uns["direct_e_regulon_metadata"]
scplus_mdata.uns["extended_e_regulon_metadata"]

#UMAP
eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,)
eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]
sc.pp.neighbors(eRegulon_gene_AUC, use_rep = "X")
sc.tl.umap(eRegulon_gene_AUC)

color_dict = {
    'AntIDIDTT': "#065143",
    'MMN': "#70B77E",
    'NPC': "#E0A890",
    'NeuralPro': "#F56476",
    'NeuralProNeurog2': "#CE1483",
    'PMN': "#053C5E" ,
    'PVHSON': "#38A3A5",
    'Prethalamus': "#80ED99",
    'SMN': "#c75c5a",
    'TuberalARCVMH': "#919191"}

sc.pl.umap(eRegulon_gene_AUC,
           color = "scRNA_counts:Cluster_Pass2", palette=color_dict)
sc.pl.umap(eRegulon_gene_AUC, color = "scRNA_counts:Cluster_Pass2")


#plot
from scenicplus.RSS import (regulon_specificity_scores, plot_rss)
rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scRNA_counts:Cluster_Pass2",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC"])
plot_rss(
    data_matrix = rss,
    top_n = 10,
    num_columns = 5)

#plot
import matplotlib.pyplot as plt
ax = sc.pl.umap(eRegulon_gene_AUC,
           color = list(set([x for xs in [rss.loc[ct].sort_values()[0:2].index for ct in rss.index] for x in xs ])))
ax.figure.savefig('umap_output.tiff', format='tiff', dpi=300)
plt.savefig('umap_output.tiff', format='tiff', dpi=300)
plt.show()

genes = ['Arx_direct_+/+_(192g)',
         'Barhl1_direct_+/+_(112g)',
         'Dlx1_direct_+/+_(296g)',
         'Foxa1_direct_+/+_(49g)',
         'Hmx2_direct_+/+_(153g)',
         'Isl1_direct_+/+_(267g)',
         'Lmx1a_direct_+/+_(118g)']
sc.pl.umap(eRegulon_gene_AUC,
           color = genes, save = True)

genes = ['Arx_direct_+/+_(192g)']
sc.pl.violin(eRegulon_gene_AUC, keys=genes,
             groupby='scRNA_counts:Cluster_Pass2')
sc.pl.violin(eRegulon_gene_AUC, keys=genes,
             groupby='scRNA_counts:Clusters', save = "test.svg")

#stacked plot
sc.pl.stacked_violin(eRegulon_gene_AUC, var_names=genes,
             groupby='scRNA_counts:Cluster_Pass2', dendrogram=False)



#plotting
# Create a new figure with specific size
fig, ax = plt.subplots(figsize=(10, 6))  # Adjust the size as needed
# Generate the plot and get the axis object
sc.pl.violin(eRegulon_gene_AUC, keys='Arx_direct_+/+_(170g)',
             groupby='scRNA_counts:Clusters', ax=ax, show=False)
# Rotate the x-axis labels
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
# Adjust layout to ensure labels are not cut off
plt.tight_layout()
# Save the plot with adjusted layout
fig.savefig("test.svg", bbox_inches='tight')
# Show the plot
plt.show()
# Close the figure to release resources
plt.close(fig)



#dotplot
from scenicplus.plotting.dotplot import heatmap_dotplot
heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scRNA_counts:Cluster_Pass2",
    group_variable_order=['AntIDIDTT', 'MMN', 'NPC', 'NeuralPro',
                          'NeuralProNeurog2', 'PMN', 'PVHSON', 'Prethalamus', 'SMN', 'TuberalARCVMH'],
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    scale_size_matrix=True,
    scale_color_matrix=True,
    orientation = "horizontal",
    figsize = (24, 10),  # Reduced to be within the default limit
    save = "heatmap_plot.png")


from scenicplus.plotting.dotplot import heatmap_dotplot
heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scRNA_counts:Clusters",
    group_variable_order=['ControlAntIDIDTT', 'MutantAntIDIDTT',
                          'ControlMMN', 'MutantMMN',
                          'ControlNPC', 'MutantNPC',
                          'ControlNeuralPro', 'MutantNeuralPro',
           'ControlNeuralProNeurog2', 'MutantNeuralProNeurog2',
           'ControlPMN', 'MutantPMN',
           'ControlPVHSON', 'MutantPVHSON',
           'ControlPrethalamus', 'MutantPrethalamus',
           'ControlSMN', 'MutantSMN',
           'ControlTuberalARCVMH', 'MutantTuberalARCVMH'],
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    scale_size_matrix=True,
    scale_color_matrix=True,
    orientation = "horizontal",
    figsize = (24, 10),  # Reduced to be within the default limit
    save = "heatmap_plot.eps")

