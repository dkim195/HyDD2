import os
os.chdir("/media/thomaskim/Data/Xenium")
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

astrcyte = ["Cdh20","Gfap","Gli3","Id2","Ntsr2"]
astrcyte_tanycyte = ["Aqp4","Glul"]
endothelial = ["Cd93","Cldn5","Ly6a","Car4","Fn1"]
tanycyte = ["Col23a1","Frzb","Rax","Slc17a8"]
microglia = ["Cd300c2","Cd68","Trem2","Siglech"]
oligo =  ["Opalin","Sox10","Gjc3","Gpr17","Sema3d"]
pericytes = ["Acta2","Cspg4", "Pln", "Ano1"]
VLMC = ["Aldh1a2","Col1a1","Dcn", "Pdgfra"]
sc.pl.violin(adata, keys = astrcyte, groupby="leiden")
sc.pl.violin(adata, keys = astrcyte_tanycyte, groupby="leiden")
sc.pl.violin(adata, keys = endothelial , groupby="leiden")
sc.pl.violin(adata, keys = tanycyte, groupby="leiden")
sc.pl.violin(adata, keys = microglia, groupby="leiden")
sc.pl.violin(adata, keys = oligo, groupby="leiden")
sc.pl.violin(adata, keys = pericytes, groupby="leiden")
sc.pl.violin(adata, keys = VLMC, groupby="leiden")

genes = ["Rax"]
sq.pl.spatial_scatter(
    adata, library_id="spatial",
    color=genes, shape=None,
    size=2, cmap="Reds", img=False, figsize=(12, 8))

dedf[dedf['group'] == '3']
sq.pl.spatial_scatter(adata,
                      shape=None, groups = ["3"],
                     color="leiden", library_id="spatial")

import os

def analyze_and_plot_groups(adata, dedf, start_group, end_group):
    # Ensure the output directory exists
    output_dir = 'png'
    os.makedirs(output_dir, exist_ok=True)
    
    # Loop through each group
    for group_id in map(str, range(start_group, end_group + 1)):
        # Filter the DataFrame for the specified group
        filtered_df = dedf[dedf['group'] == group_id]
        print(f"Data for group {group_id}:")
        print(filtered_df)
        
        # Plotting the spatial scatter for the group
        print(f"Plotting spatial scatter for group {group_id}...")
        plot_path = os.path.join(output_dir, f'group_{group_id}_scatter.png')
        sq.pl.spatial_scatter(adata, shape=None, groups=[group_id], color="leiden", library_id="spatial", save=plot_path)
        print(f"Plot saved to {plot_path}")

# Usage
analyze_and_plot_groups(adata, dedf, 0, 22)


#HET
#Region1
adata = sc.read_h5ad(filename="OUTPUT/P21_HET_REGION1.h5ad")
sq.pl.spatial_scatter(adata,
                      shape=None, color="leiden", library_id="spatial")

sc.tl.rank_genes_groups(adata, groupby="leiden", use_raw=False)
group_names = list(map(str, range(23))) 
dedf = sc.get.rank_genes_groups_df(adata, group = group_names)

def analyze_and_plot_group(adata, dedf, group_id):
    # Filter the DataFrame for the specified group
    filtered_df = dedf[dedf['group'] == group_id]
    print(filtered_df)
    sq.pl.spatial_scatter(adata, shape=None, groups=[group_id], color="leiden", library_id="spatial")
    
# Specify the group ID as a string if the group column is string-typed
group_id = '0'  
analyze_and_plot_group(adata, dedf, group_id)

old_to_new = {
"0":'Excitatory_Telencephalon_1',
"1":'Oligodendrocyte',
"2":'Endothelial',
"3":'Inhibitory_1',
"4":'Excitatory_Diencephalon',
"5":'Astrocyte_1',
"6":'Astrocyte_2',
"7":'Astrocyte_3',
"8":'Endothelial',
"9":'OPC_1',
"10":'Hippocampus',
"11":'Microglia',
"12":'Excitatory_Telencephalon_2',
"13":'Inhibitory_Telencephalon_1',
"14":'Hippocampus',
"15":'Ependymal',
"16":'Thalamus_1',
"17":'Thalamus_2',
"18":'OPC_2',
"19":'Inhibitory_2',
"20":'VSMC_Pericyte',
"21":'Inhibitory_Telencephalon_2',
"22":'NTS_Neuron'}
adata.obs['annotation'] = (
adata.obs['leiden']
.map(old_to_new)
.astype('category')
)
if 'annotation_colors' in adata.uns:
    adata.uns.pop('annotation_colors')
    
annotation_colors = {
    'Excitatory_Telencephalon_1': '#e6194B',
    'Oligodendrocyte': '#3cb44b',
    'Endothelial': '#ffe119',
    'Inhibitory_1': '#4363d8',
    'Excitatory_Diencephalon': '#f58231',
    'Astrocyte_1': '#911eb4',
    'Astrocyte_2': '#42d4f4',
    'Astrocyte_3': '#f032e6',
    'OPC_1': '#bfef45',
    'Hippocampus': '#fabebe',
    'Microglia': '#469990',
    'Excitatory_Telencephalon_2': '#e6beff',
    'Inhibitory_Telencephalon_1': '#9A6324',
    'Ependymal': '#fffac8',
    'Thalamus_1': '#800000',
    'Thalamus_2': '#aaffc3',
    'OPC_2': '#808000',
    'Inhibitory_2': '#ffd8b1',
    'VSMC_Pericyte': '#000075',
    'Inhibitory_Telencephalon_2': '#a9a9a9', 
    'NTS_Neuron': '#000000'
}

import matplotlib.colors as mcolors
categories = adata.obs['annotation'].cat.categories
palette_list = [annotation_colors[category] for category in categories if category in annotation_colors]
listed_colormap = mcolors.ListedColormap(palette_list)

sq.pl.spatial_scatter(adata,
                      shape=None, color="annotation",
                      library_id="spatial", legend_loc="right",
                      legend_fontsize = "xx-small",
                      palette = listed_colormap, figsize=(20, 10))

sq.gr.nhood_enrichment(adata, cluster_key="annotation")
sq.pl.nhood_enrichment(adata, cluster_key="annotation", figsize=(5, 5))

adata.write("OUTPUT/P21_HET_REGION1.h5ad")







#REGION2
adata = sc.read_h5ad(filename="OUTPUT/P21_HET_REGION2.h5ad")
sq.pl.spatial_scatter(adata,
                      shape=None, color="leiden", library_id="spatial")

sc.tl.rank_genes_groups(adata, groupby="leiden", use_raw=False)
group_names = list(map(str, range(25))) 
dedf = sc.get.rank_genes_groups_df(adata, group = group_names)

def analyze_and_plot_group(adata, dedf, group_id):
    # Filter the DataFrame for the specified group
    filtered_df = dedf[dedf['group'] == group_id]
    print(filtered_df)
    sq.pl.spatial_scatter(adata, shape=None, groups=[group_id], color="leiden", library_id="spatial")
    
# Specify the group ID as a string if the group column is string-typed
group_id = '26'  
analyze_and_plot_group(adata, dedf, group_id)

old_to_new = {
"0":'Oligodendrocyte',
"1":'Excitatory_Telencephalon_1',
"2":'Astrocyte_3',
"3":'Endothelial',
"4":'Inhibitory_1',
"5":'Excitatory_Telencephalon_3',
"6":'Excitatory_Telencephalon_2',
"7":'Excitatory_Diencephalon',
"8":'Endothelial',
"9":'OPC_1',
"10":'Hippocampus',
"11":'Microglia',
"12":'Astrocyte_1',
"13":'Astrocyte_1',
"14":'Hippocampus',
"15":'Thalamus_1',
"16":'Inhibitory_2',
"17":'Inhibitory_Telencephalon_1',
"18":'Ependymal',
"19":'OPC_2',
"20":'Excitatory_Diencephalon',
"21":'Endothelial',
"22":'Inhibitory_Telencephalon_3',
"23":'Thalamus_3',
"24":'UNKNOWN',
"25":'UNKNOWN',
"26":'UNKNOWN'}
adata.obs['annotation'] = (
adata.obs['leiden']
.map(old_to_new)
.astype('category')
)
if 'annotation_colors' in adata.uns:
    adata.uns.pop('annotation_colors')
    
annotation_colors = {
    'Excitatory_Telencephalon_1': '#e6194B',
    'Oligodendrocyte': '#3cb44b',
    'Endothelial': '#ffe119',
    'Inhibitory_1': '#4363d8',
    'Excitatory_Diencephalon': '#f58231',
    'Astrocyte_1': '#911eb4',
    'Astrocyte_2': '#42d4f4',
    'Astrocyte_3': '#f032e6',
    'OPC_1': '#bfef45',
    'Hippocampus': '#fabebe',
    'Microglia': '#469990',
    'Excitatory_Telencephalon_2': '#e6beff',
    'Inhibitory_Telencephalon_1': '#9A6324',
    'Ependymal': '#fffac8',
    'Thalamus_1': '#800000',
    'Thalamus_2': '#aaffc3',
    'OPC_2': '#808000',
    'Inhibitory_2': '#ffd8b1',
    'VSMC_Pericyte': '#000075',
    'Inhibitory_Telencephalon_2': '#a9a9a9', 
    'NTS_Neuron': '#000000',
    'Excitatory_Telencephalon_3': '#40e0d0',
    'Inhibitory_Telencephalon_3': '#DC143C',
    'Thalamus_3':'#00BFFF', 
    'UNKNOWN': '#D3D3D3'
}

import matplotlib.colors as mcolors
categories = adata.obs['annotation'].cat.categories
palette_list = [annotation_colors[category] for category in categories if category in annotation_colors]
listed_colormap = mcolors.ListedColormap(palette_list)

sq.pl.spatial_scatter(adata,
                      shape=None, color="annotation",
                      library_id="spatial", legend_loc="right",
                      legend_fontsize = "xx-small",
                      palette = listed_colormap, figsize=(20, 10))

sq.gr.nhood_enrichment(adata, cluster_key="annotation")
sq.pl.nhood_enrichment(adata, cluster_key="annotation", figsize=(5, 5))

adata.write("OUTPUT/P21_HET_REGION2.h5ad")





#HET
#REGION3
adata = sc.read_h5ad(filename="OUTPUT/P21_HET_REGION3.h5ad")
sq.pl.spatial_scatter(adata,
                      shape=None, color="leiden", library_id="spatial")

sc.tl.rank_genes_groups(adata, groupby="leiden", use_raw=False)
group_names = list(map(str, range(26))) 
dedf = sc.get.rank_genes_groups_df(adata, group = group_names)

def analyze_and_plot_group(adata, dedf, group_id):
    # Filter the DataFrame for the specified group
    filtered_df = dedf[dedf['group'] == group_id]
    print(filtered_df)
    sq.pl.spatial_scatter(adata, shape=None, groups=[group_id], color="leiden", library_id="spatial")
    
# Specify the group ID as a string if the group column is string-typed
group_id = '26'  
analyze_and_plot_group(adata, dedf, group_id)

old_to_new = {
"0":'Excitatory_Telencephalon_1',
"1":'Inhibitory_1',
"2":'Endothelial',
"3":'Oligodendrocyte',
"4":'Astrocyte_2',
"5":'Endothelial',
"6":'OPC_1',
"7":'Microglia',
"8":'Oligodendrocyte',
"9":'Excitatory_Diencephalon',
"10":'Inhibitory_Telencephalon_1',
"11":'Astrocyte_1',
"12":'Excitatory_Telencephalon_2',
"13":'Hippocampus',
"14":'Excitatory_Telencephalon_3',
"15":'Thalamus_1',
"16":'Thalamus_3',
"17":'Hippocampus',
"18":'Astrocyte_1',
"19":'Inhibitory_2',
"20":'OPC_2',
"21":'Ependymal',
"22":'Excitatory_Telencephalon_2',
"23":'UNKNOWN',
"24":'Endothelial',
"25":'UNKNOWN',
"26":'UNKNOWN'}
adata.obs['annotation'] = (
adata.obs['leiden']
.map(old_to_new)
.astype('category')
)
if 'annotation_colors' in adata.uns:
    adata.uns.pop('annotation_colors')
    
annotation_colors = {
    'Excitatory_Telencephalon_1': '#e6194B',
    'Oligodendrocyte': '#3cb44b',
    'Endothelial': '#ffe119',
    'Inhibitory_1': '#4363d8',
    'Excitatory_Diencephalon': '#f58231',
    'Astrocyte_1': '#911eb4',
    'Astrocyte_2': '#42d4f4',
    'Astrocyte_3': '#f032e6',
    'OPC_1': '#bfef45',
    'Hippocampus': '#fabebe',
    'Microglia': '#469990',
    'Excitatory_Telencephalon_2': '#e6beff',
    'Inhibitory_Telencephalon_1': '#9A6324',
    'Ependymal': '#fffac8',
    'Thalamus_1': '#800000',
    'Thalamus_2': '#aaffc3',
    'OPC_2': '#808000',
    'Inhibitory_2': '#ffd8b1',
    'VSMC_Pericyte': '#000075',
    'Inhibitory_Telencephalon_2': '#a9a9a9', 
    'NTS_Neuron': '#000000',
    'Excitatory_Telencephalon_3': '#40e0d0',
    'Inhibitory_Telencephalon_3': '#DC143C',
    'Thalamus_3':'#00BFFF', 
    'UNKNOWN': '#D3D3D3'
}

import matplotlib.colors as mcolors
categories = adata.obs['annotation'].cat.categories
palette_list = [annotation_colors[category] for category in categories if category in annotation_colors]
listed_colormap = mcolors.ListedColormap(palette_list)

sq.pl.spatial_scatter(adata,
                      shape=None, color="annotation",
                      library_id="spatial", legend_loc="right",
                      legend_fontsize = "xx-small",
                      palette = listed_colormap, figsize=(20, 10))

sq.gr.nhood_enrichment(adata, cluster_key="annotation")
sq.pl.nhood_enrichment(adata, cluster_key="annotation", figsize=(5, 5))

adata.write("OUTPUT/P21_HET_REGION3.h5ad")



#HET
#REGION4
adata = sc.read_h5ad(filename="OUTPUT/P21_HET_REGION4.h5ad")
sq.pl.spatial_scatter(adata,
                      shape=None, color="leiden", library_id="spatial")

sc.tl.rank_genes_groups(adata, groupby="leiden", use_raw=False)
group_names = list(map(str, range(26))) 
dedf = sc.get.rank_genes_groups_df(adata, group = group_names)

def analyze_and_plot_group(adata, dedf, group_id):
    # Filter the DataFrame for the specified group
    filtered_df = dedf[dedf['group'] == group_id]
    print(filtered_df)
    sq.pl.spatial_scatter(adata, shape=None, groups=[group_id], color="leiden", library_id="spatial")
    
# Specify the group ID as a string if the group column is string-typed
group_id = '26'  
analyze_and_plot_group(adata, dedf, group_id)

old_to_new = {
"0":'Oligodendrocyte',
"1":'Endothelial',
"2":'Inhibitory_1',
"3":'Excitatory_Telencephalon_1',
"4":'Astrocyte_2',
"5":'Inhibitory_Telencephalon_1',
"6":'Thalamus_3',
"7":'Excitatory_Telencephalon_2',
"8":'Endothelial',
"9":'OPC_1',
"10":'Microglia',
"11":'Excitatory_Telencephalon_2',
"12":'Excitatory_Diencephalon',
"13":'Inhibitory_2',
"14":'Excitatory_Telencephalon_2',
"15":'Hippocampus',
"16":'Astrocyte_1',
"17":'OPC_2',
"18":'Hippocampus',
"19":'Hippocampus',
"20":'Ependymal',
"21":'UNKNOWN',
"22":'Excitatory_Telencephalon_2',
"23":'SCN',
"24":'Endothelial',
"25":'Astrocyte_2'}
adata.obs['annotation'] = (
adata.obs['leiden']
.map(old_to_new)
.astype('category')
)
if 'annotation_colors' in adata.uns:
    adata.uns.pop('annotation_colors')
    
annotation_colors = {
    'Excitatory_Telencephalon_1': '#e6194B',
    'Oligodendrocyte': '#3cb44b',
    'Endothelial': '#ffe119',
    'Inhibitory_1': '#4363d8',
    'Excitatory_Diencephalon': '#f58231',
    'Astrocyte_1': '#911eb4',
    'Astrocyte_2': '#42d4f4',
    'Astrocyte_3': '#f032e6',
    'OPC_1': '#bfef45',
    'Hippocampus': '#fabebe',
    'Microglia': '#469990',
    'Excitatory_Telencephalon_2': '#e6beff',
    'Inhibitory_Telencephalon_1': '#9A6324',
    'Ependymal': '#fffac8',
    'Thalamus_1': '#800000',
    'Thalamus_2': '#aaffc3',
    'OPC_2': '#808000',
    'Inhibitory_2': '#ffd8b1',
    'VSMC_Pericyte': '#000075',
    'Inhibitory_Telencephalon_2': '#a9a9a9', 
    'NTS_Neuron': '#000000',
    'Excitatory_Telencephalon_3': '#40e0d0',
    'Inhibitory_Telencephalon_3': '#DC143C',
    'Thalamus_3':'#00BFFF', 
    'UNKNOWN': '#D3D3D3',
    'SCN': '#00FF00'
}

import matplotlib.colors as mcolors
categories = adata.obs['annotation'].cat.categories
palette_list = [annotation_colors[category] for category in categories if category in annotation_colors]
listed_colormap = mcolors.ListedColormap(palette_list)

sq.pl.spatial_scatter(adata,
                      shape=None, color="annotation",
                      library_id="spatial", legend_loc="right",
                      legend_fontsize = "xx-small",
                      palette = listed_colormap, figsize=(20, 10))

sq.gr.nhood_enrichment(adata, cluster_key="annotation")
sq.pl.nhood_enrichment(adata, cluster_key="annotation", figsize=(5, 5))

adata.write("OUTPUT/P21_HET_REGION4.h5ad")





#HOMO
#REGION1
adata = sc.read_h5ad(filename="OUTPUT/P21_HOMO_REGION1.h5ad")
sq.pl.spatial_scatter(adata,
                      shape=None, color="leiden", library_id="spatial")

sc.tl.rank_genes_groups(adata, groupby="leiden", use_raw=False)
group_names = list(map(str, range(23))) 
dedf = sc.get.rank_genes_groups_df(adata, group = group_names)

def analyze_and_plot_group(adata, dedf, group_id):
    # Filter the DataFrame for the specified group
    filtered_df = dedf[dedf['group'] == group_id]
    print(filtered_df)
    sq.pl.spatial_scatter(adata, shape=None, groups=[group_id], color="leiden", library_id="spatial")
    
# Specify the group ID as a string if the group column is string-typed
group_id = '24'  
analyze_and_plot_group(adata, dedf, group_id)

old_to_new = {
"0":'Astrocyte_2',
"1":'Endothelial',
"2":'Oligodendrocyte',
"3":'Excitatory_Telencephalon_1',
"4":'Excitatory_Diencephalon',
"5":'Endothelial',
"6":'Excitatory_Telencephalon_2',
"7":'Excitatory_Telencephalon_3',
"8":'Inhibitory_Telencephalon_3',
"9":'OPC_1',
"10":'Microglia',
"11":'Astrocyte_3',
"12":'Inhibitory_1',
"13":'Excitatory_Telencephalon_2',
"14":'Hippocampus',
"15":'Thalamus_1',
"16":'Hippocampus',
"17":'Ependymal',
"18":'OPC_2',
"19":'Endothelial',
"20":'Inhibitory_Diencephalon_1',
"21":'Excitatory_Telencephalon_2',
"22":'UNKNOWN',
"23":'UNKNOWN',
"24":'UNKNOWN'}
adata.obs['annotation'] = (
adata.obs['leiden']
.map(old_to_new)
.astype('category')
)
if 'annotation_colors' in adata.uns:
    adata.uns.pop('annotation_colors')
    
annotation_colors = {
    'Excitatory_Telencephalon_1': '#e6194B',
    'Oligodendrocyte': '#3cb44b',
    'Endothelial': '#ffe119',
    'Inhibitory_1': '#4363d8',
    'Excitatory_Diencephalon': '#f58231',
    'Astrocyte_1': '#911eb4',
    'Astrocyte_2': '#42d4f4',
    'Astrocyte_3': '#f032e6',
    'OPC_1': '#bfef45',
    'Hippocampus': '#fabebe',
    'Microglia': '#469990',
    'Excitatory_Telencephalon_2': '#e6beff',
    'Inhibitory_Telencephalon_1': '#9A6324',
    'Ependymal': '#fffac8',
    'Thalamus_1': '#800000',
    'Thalamus_2': '#aaffc3',
    'OPC_2': '#808000',
    'Inhibitory_2': '#ffd8b1',
    'VSMC_Pericyte': '#000075',
    'Inhibitory_Telencephalon_2': '#a9a9a9', 
    'NTS_Neuron': '#000000',
    'Excitatory_Telencephalon_3': '#40e0d0',
    'Inhibitory_Telencephalon_3': '#DC143C',
    'Thalamus_3':'#00BFFF', 
    'UNKNOWN': '#D3D3D3',
    'Inhibitory_Diencephalon_1': '#00FF00'
}

import matplotlib.colors as mcolors
categories = adata.obs['annotation'].cat.categories
palette_list = [annotation_colors[category] for category in categories if category in annotation_colors]
listed_colormap = mcolors.ListedColormap(palette_list)

sq.pl.spatial_scatter(adata,
                      shape=None, color="annotation",
                      library_id="spatial", legend_loc="right",
                      legend_fontsize = "xx-small",
                      palette = listed_colormap, figsize=(20, 10))

sq.gr.nhood_enrichment(adata, cluster_key="annotation")
sq.pl.nhood_enrichment(adata, cluster_key="annotation", figsize=(5, 5))

adata.write("OUTPUT/P21_HOMO_REGION1.h5ad")




#REGION2
adata = sc.read_h5ad(filename="OUTPUT/P21_HOMO_REGION2.h5ad")
sq.pl.spatial_scatter(adata,
                      shape=None, color="leiden", library_id="spatial")

sc.tl.rank_genes_groups(adata, groupby="leiden", use_raw=False)
group_names = list(map(str, range(26))) 
dedf = sc.get.rank_genes_groups_df(adata, group = group_names)

def analyze_and_plot_group(adata, dedf, group_id):
    # Filter the DataFrame for the specified group
    filtered_df = dedf[dedf['group'] == group_id]
    print(filtered_df)
    sq.pl.spatial_scatter(adata, shape=None, groups=[group_id], color="leiden", library_id="spatial")
    
# Specify the group ID as a string if the group column is string-typed
group_id = '27'  
analyze_and_plot_group(adata, dedf, group_id)

old_to_new = {
"0":'Endothelial',
"1":'Oligodendrocyte',
"2":'Astrocyte_3',
"3":'Excitatory_Telencephalon_1',
"4":'Excitatory_Telencephalon_3',
"5":'Inhibitory_1',
"6":'Thalamus_3',
"7":'Endothelial',
"8":'Astrocyte_1',
"9":'OPC_1',
"10":'Hippocampus',
"11":'Excitatory_Telencephalon_2',
"12":'Microglia',
"13":'Inhibitory_Diencephalon_1',
"14":'Endothelial',
"15":'Thalamus_1',
"16":'OPC_2',
"17":'Hippocampus',
"18":'Inhibitory_2',
"19":'Thalamus_2',
"20":'Ependymal',
"21":'Excitatory_Telencephalon_2',
"22":'Inhibitory_Telencephalon_1',
"23":'UNKNOWN',
"24":'UNKNOWN',
"25":'MUTANT_UNKNOWN',
"26":'UNKNOWN',
"26":'UNKNOWN',
"27":'UNKNOWN'}
adata.obs['annotation'] = (
adata.obs['leiden']
.map(old_to_new)
.astype('category')
)
if 'annotation_colors' in adata.uns:
    adata.uns.pop('annotation_colors')
    
annotation_colors = {
    'Excitatory_Telencephalon_1': '#e6194B',
    'Oligodendrocyte': '#3cb44b',
    'Endothelial': '#ffe119',
    'Inhibitory_1': '#4363d8',
    'Excitatory_Diencephalon': '#f58231',
    'Astrocyte_1': '#911eb4',
    'Astrocyte_2': '#42d4f4',
    'Astrocyte_3': '#f032e6',
    'OPC_1': '#bfef45',
    'Hippocampus': '#fabebe',
    'Microglia': '#469990',
    'Excitatory_Telencephalon_2': '#e6beff',
    'Inhibitory_Telencephalon_1': '#9A6324',
    'Ependymal': '#fffac8',
    'Thalamus_1': '#800000',
    'Thalamus_2': '#aaffc3',
    'OPC_2': '#808000',
    'Inhibitory_2': '#ffd8b1',
    'VSMC_Pericyte': '#000075',
    'Inhibitory_Telencephalon_2': '#a9a9a9', 
    'NTS_Neuron': '#000000',
    'Excitatory_Telencephalon_3': '#40e0d0',
    'Inhibitory_Telencephalon_3': '#DC143C',
    'Thalamus_3':'#00BFFF', 
    'UNKNOWN': '#D3D3D3',
    'Inhibitory_Diencephalon_1': '#00FF00',
    'MUTANT_UNKNOWN':'#FF2400'
}

import matplotlib.colors as mcolors
categories = adata.obs['annotation'].cat.categories
palette_list = [annotation_colors[category] for category in categories if category in annotation_colors]
listed_colormap = mcolors.ListedColormap(palette_list)

sq.pl.spatial_scatter(adata,
                      shape=None, color="annotation",
                      library_id="spatial", legend_loc="right",
                      legend_fontsize = "xx-small",
                      palette = listed_colormap, figsize=(20, 10))

sq.gr.nhood_enrichment(adata, cluster_key="annotation")
sq.pl.nhood_enrichment(adata, cluster_key="annotation", figsize=(5, 5))

adata.write("OUTPUT/P21_HOMO_REGION2.h5ad")





#REGION3
adata = sc.read_h5ad(filename="OUTPUT/P21_HOMO_REGION3.h5ad")
sq.pl.spatial_scatter(adata,
                      shape=None, color="leiden", library_id="spatial")

sc.tl.rank_genes_groups(adata, groupby="leiden", use_raw=False)
group_names = list(map(str, range(26))) 
dedf = sc.get.rank_genes_groups_df(adata, group = group_names)

def analyze_and_plot_group(adata, dedf, group_id):
    # Filter the DataFrame for the specified group
    filtered_df = dedf[dedf['group'] == group_id]
    print(filtered_df)
    sq.pl.spatial_scatter(adata, shape=None, groups=[group_id], color="leiden", library_id="spatial")
    
# Specify the group ID as a string if the group column is string-typed
group_id = '27'  
analyze_and_plot_group(adata, dedf, group_id)

old_to_new = {
"0":'Oligodendrocyte',
"1":'Astrocyte_3',
"2":'Endothelial',
"3":'Excitatory_Telencephalon_1',
"4":'Inhibitory_1',
"5":'Endothelial',
"6":'Astrocyte_1',
"7":'OPC_1',
"8":'Inhibitory_Diencephalon_1',
"9":'Microglia',
"10":'Excitatory_Telencephalon_2',
"11":'Thalamus_3',
"12":'Thalamus_2',
"13":'Excitatory_Telencephalon_3',
"14":'Excitatory_Telencephalon_2',
"15":'OPC_2',
"16":'Ependymal',
"17":'Inhibitory_Telencephalon_1',
"18":'Hippocampus',
"19":'Endothelial',
"20":'Hippocampus',
"21":'Inhibitory_Telencephalon_3',
"22":'Excitatory_Telencephalon_2',
"23":'UNKNOWN',
"24":'MUTANT_UNKNOWN',
"25":'UNKNOWN',
"26":'UNKNOWN',
"26":'UNKNOWN',
"27":'UNKNOWN'}
adata.obs['annotation'] = (
adata.obs['leiden']
.map(old_to_new)
.astype('category')
)
if 'annotation_colors' in adata.uns:
    adata.uns.pop('annotation_colors')
    
annotation_colors = {
    'Excitatory_Telencephalon_1': '#e6194B',
    'Oligodendrocyte': '#3cb44b',
    'Endothelial': '#ffe119',
    'Inhibitory_1': '#4363d8',
    'Excitatory_Diencephalon': '#f58231',
    'Astrocyte_1': '#911eb4',
    'Astrocyte_2': '#42d4f4',
    'Astrocyte_3': '#f032e6',
    'OPC_1': '#bfef45',
    'Hippocampus': '#fabebe',
    'Microglia': '#469990',
    'Excitatory_Telencephalon_2': '#e6beff',
    'Inhibitory_Telencephalon_1': '#9A6324',
    'Ependymal': '#fffac8',
    'Thalamus_1': '#800000',
    'Thalamus_2': '#aaffc3',
    'OPC_2': '#808000',
    'Inhibitory_2': '#ffd8b1',
    'VSMC_Pericyte': '#000075',
    'Inhibitory_Telencephalon_2': '#a9a9a9', 
    'NTS_Neuron': '#000000',
    'Excitatory_Telencephalon_3': '#40e0d0',
    'Inhibitory_Telencephalon_3': '#DC143C',
    'Thalamus_3':'#00BFFF', 
    'UNKNOWN': '#D3D3D3',
    'Inhibitory_Diencephalon_1': '#00FF00',
    'MUTANT_UNKNOWN':'#FF2400'
}

import matplotlib.colors as mcolors
categories = adata.obs['annotation'].cat.categories
palette_list = [annotation_colors[category] for category in categories if category in annotation_colors]
listed_colormap = mcolors.ListedColormap(palette_list)

sq.pl.spatial_scatter(adata,
                      shape=None, color="annotation",
                      library_id="spatial", legend_loc="right",
                      legend_fontsize = "xx-small",
                      palette = listed_colormap, figsize=(20, 10))

sq.gr.nhood_enrichment(adata, cluster_key="annotation")
sq.pl.nhood_enrichment(adata, cluster_key="annotation", figsize=(5, 5))

adata.write("OUTPUT/P21_HOMO_REGION3.h5ad")



#REGION4
adata = sc.read_h5ad(filename="OUTPUT/P21_HOMO_REGION4.h5ad")
sq.pl.spatial_scatter(adata,
                      shape=None, color="leiden", library_id="spatial")

sc.tl.rank_genes_groups(adata, groupby="leiden", use_raw=False)
group_names = list(map(str, range(24))) 
dedf = sc.get.rank_genes_groups_df(adata, group = group_names)

def analyze_and_plot_group(adata, dedf, group_id):
    # Filter the DataFrame for the specified group
    filtered_df = dedf[dedf['group'] == group_id]
    print(filtered_df)
    sq.pl.spatial_scatter(adata, shape=None, groups=[group_id], color="leiden", library_id="spatial")
    
# Specify the group ID as a string if the group column is string-typed
group_id = '25'  
analyze_and_plot_group(adata, dedf, group_id)

old_to_new = {
"0":'Astrocyte_3',
"1":'Excitatory_Telencephalon_1',
"2":'Oligodendrocyte',
"3":'Endothelial',
"4":'Inhibitory_1',
"5":'Endothelial',
"6":'OPC_1',
"7":'Microglia',
"8":'Excitatory_Diencephalon',
"9":'Inhibitory_Telencephalon_1',
"10":'Hippocampus',
"11":'Thalamus_2',
"12":'Excitatory_Telencephalon_3',
"13":'Excitatory_Telencephalon_2',
"14":'Excitatory_Telencephalon_2',
"15":'OPC_2',
"16":'Hippocampus',
"17":'Thalamus_3',
"18":'Hippocampus',
"19":'Endothelial',
"20":'Ependymal',
"21":'Inhibitory_Telencephalon_3',
"22":'UNKNOWN',
"23":'Hippocampus',
"24":'UNKNOWN',
"25":'UNKNOWN'}
adata.obs['annotation'] = (
adata.obs['leiden']
.map(old_to_new)
.astype('category')
)
if 'annotation_colors' in adata.uns:
    adata.uns.pop('annotation_colors')
    
annotation_colors = {
    'Excitatory_Telencephalon_1': '#e6194B',
    'Oligodendrocyte': '#3cb44b',
    'Endothelial': '#ffe119',
    'Inhibitory_1': '#4363d8',
    'Excitatory_Diencephalon': '#f58231',
    'Astrocyte_1': '#911eb4',
    'Astrocyte_2': '#42d4f4',
    'Astrocyte_3': '#f032e6',
    'OPC_1': '#bfef45',
    'Hippocampus': '#fabebe',
    'Microglia': '#469990',
    'Excitatory_Telencephalon_2': '#e6beff',
    'Inhibitory_Telencephalon_1': '#9A6324',
    'Ependymal': '#fffac8',
    'Thalamus_1': '#800000',
    'Thalamus_2': '#aaffc3',
    'OPC_2': '#808000',
    'Inhibitory_2': '#ffd8b1',
    'VSMC_Pericyte': '#000075',
    'Inhibitory_Telencephalon_2': '#a9a9a9', 
    'NTS_Neuron': '#000000',
    'Excitatory_Telencephalon_3': '#40e0d0',
    'Inhibitory_Telencephalon_3': '#DC143C',
    'Thalamus_3':'#00BFFF', 
    'UNKNOWN': '#D3D3D3',
    'Inhibitory_Diencephalon_1': '#00FF00',
    'MUTANT_UNKNOWN':'#FF2400'
}

import matplotlib.colors as mcolors
categories = adata.obs['annotation'].cat.categories
palette_list = [annotation_colors[category] for category in categories if category in annotation_colors]
listed_colormap = mcolors.ListedColormap(palette_list)

sq.pl.spatial_scatter(adata,
                      shape=None, color="annotation",
                      library_id="spatial", legend_loc="right",
                      legend_fontsize = "xx-small",
                      palette = listed_colormap, figsize=(20, 10))

sq.gr.nhood_enrichment(adata, cluster_key="annotation")
sq.pl.nhood_enrichment(adata, cluster_key="annotation", figsize=(5, 5))

adata.write("OUTPUT/P21_HOMO_REGION4.h5ad")
