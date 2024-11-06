#!/bin/bash
#SBATCH --account=Hypothalamus
#SBATCH -c 20
#SBATCH --mem=512g
#SBATCH --partition=normal
#SBATCH --time=2:00:00
#SBATCH --mail-type=END
#SBATCH --job-name=DLX_Ctrl_scenicplus_pipeline
#SBATCH --mail-user=tkim@dandrite.au.dk

# Print the start time of the job
cd /faststorage/project/Hypothalamus/SCENIC/
echo "Job started at: $(date)" >> DLX_Ctrl_scenicplus_pipeline.log


# Proxy configuration
export http_proxy="http://proxy-default:3128"
export https_proxy="http://proxy-default:3128"
echo "Proxy configured for HTTP and HTTPS" >> DLX_Ctrl_scenicplus_pipeline.log

# Initialize Conda
echo "Initializing Conda" >> DLX_Ctrl_scenicplus_pipeline.log
source /home/tkim/miniforge3/etc/profile.d/conda.sh

# Activate the Conda environment
echo "Activating Conda environment scenicplus" >> DLX_Ctrl_scenicplus_pipeline.log
conda activate scenicplus

# Define NCBI IDs
NCBI_GENOME_ID="52"  # Mus musculus example
NCBI_ASSEMBLY_ID="327618"  # Specific assembly for GRCm38.p6
echo "Using hardcoded NCBI Genome ID: $NCBI_GENOME_ID" >> DLX_Ctrl_scenicplus_pipeline.log
echo "Using hardcoded NCBI Assembly ID: $NCBI_ASSEMBLY_ID" >> DLX_Ctrl_scenicplus_pipeline.log

#Update dask to the required version
#echo "Updating dask to version 2024.05" >> DLX_Ctrl_scenicplus_pipeline.log
#pip uninstall -y dask  # Added -y to avoid confirmation prompt
#pip install dask==2024.2.1

# Check if the conda command is available
if ! command -v conda &> /dev/null; then
    echo "Conda command could not be found" >> DLX_Ctrl_scenicplus_pipeline.log
    exit 1
fi

# Check if the environment was activated successfully
if [ $? -ne 0 ]; then
    echo "Failed to activate Conda environment scenicplus" >> DLX_Ctrl_scenicplus_pipeline.log
    exit 1
fi

# Log the currently active Conda environment
echo "Current Conda environment:" >> DLX_Ctrl_scenicplus_pipeline.log
conda info --envs >> DLX_Ctrl_scenicplus_pipeline.log

# Change back to the working directory
echo "Changing to working directory" >> DLX_Ctrl_scenicplus_pipeline.log
cd /faststorage/project/Hypothalamus/SCENIC/

# More commands as per your original script, incorporating proxy-safe network calls
# Example: using curl with proxy settings directly
curl --proxy "http://proxy-default:3128" --verbose https://ftp.ncbi.nlm.nih.gov/ --output test_connection.log

# Continue with the rest of your script

# Process adata object with Scanpy
echo "Processing adata object with Scanpy" >> DLX_Ctrl_scenicplus_pipeline.log
python - <<EOF 2>> DLX_Ctrl_scenicplus_pipeline.log
import scanpy as sc
import logging
import pandas as pd

# Configure logging
logging.basicConfig(filename='DLX_Ctrl_scenicplus_pipeline.log', level=logging.DEBUG)

try:
    logging.info("Loading AnnData object")
    adata = sc.read_h5ad("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT/DLX/CTRL_DLX.h5ad")

    # Clean up and consolidate cluster names

    # Adjusting specific clusters in the 'Cluster_Pass2' column

    adata.obs['Cluster_Pass2'] = (
        adata.obs['Cluster_Pass2']
        .str.replace(r'\(', '', regex=True)
        .str.replace(r'\)', '', regex=True)
        .str.replace(r' ', '_', regex=True)
        .str.replace(r',', '', regex=True)
    )

    adata.obs['Cluster_Pass2'] = (
        adata.obs['Cluster_Pass2']
        .str.replace(r'Prethalamus_.+', 'Prethalamus', regex=True)  # Consolidate all 'Prethalamus_' types
        .str.replace(r'PVH_SON_.+', 'PVH_SON', regex=True)          # Consolidate all 'PVH_SON_' types
    )

    adata.obs['Cluster_Pass2'] = (
        adata.obs['Cluster_Pass2']
        .str.replace(r'_', '', regex=True) 
        .str.replace(r'-', '', regex=True)
    )

    # Store the raw counts before normalization
    logging.info("Storing raw counts")
    new_adata = adata.raw.to_adata()

    # Set 'var_names' to the index of the 'raw.var' DataFrame, which contains your gene names
    new_adata.var_names = adata.raw.var.index
    new_adata.var_names.name = None
    new_adata.raw = new_adata

    # Normalize and log transform the data
    logging.info("Normalizing and log transforming the data")
    sc.pp.normalize_total(new_adata, target_sum=1e4)
    sc.pp.log1p(new_adata)

    # Save the processed AnnData object
    logging.info("Saving the processed AnnData object")
    new_adata.write("/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT/DLX/CTRL_DLX2.h5ad")
    logging.info("Processing completed successfully")

    # Getting a new summary of the groups in the 'Cluster_Pass2' column
    new_summary = adata.obs['Cluster_Pass2'].value_counts()

    # Printing the new summary
    print(new_summary)

except Exception as e:
    logging.error(f"An error occurred: {e}")
    raise
EOF


# Check if the Python script exited successfully
if [ $? -ne 0 ]; then
    echo "Python script for processing adata failed" >> DLX_Ctrl_scenicplus_pipeline.log
    exit 1
fi

echo "Python script for processing adata completed successfully" >> DLX_Ctrl_scenicplus_pipeline.log

# Example of using curl with proxy settings directly
echo "Testing connection to NCBI." >> DLX_Ctrl_scenicplus_pipeline.log
curl --proxy "http://proxy-default:3128" --verbose -f -L https://ftp.ncbi.nlm.nih.gov/ --output test_connection.log 2>> DLX_Ctrl_scenicplus_pipeline.log
if [ $? -ne 0 ]; then
    echo "Failed to download data from NCBI" >> DLX_Ctrl_scenicplus_pipeline.log
    exit 1
fi

# Initialize Snakemake pipeline
echo "Initializing Snakemake pipeline" >> DLX_Ctrl_scenicplus_pipeline.log

# Check if Snakemake directory already exists
if [ -d "scplus_pipeline/Snakemake" ]; then
    echo "Removing existing Snakemake directory" >> DLX_Ctrl_scenicplus_pipeline.log
    rm -rf scplus_pipeline/Snakemake
fi

# Create the Snakemake directory
mkdir -p scplus_pipeline
scenicplus init_snakemake --out_dir scplus_pipeline

# Verify the Snakemake pipeline was initialized
if [ ! -d "scplus_pipeline/Snakemake" ]; then
    echo "Snakemake pipeline initialization failed" >> DLX_Ctrl_scenicplus_pipeline.log
    exit 1
fi

# Create output and temporary directories
echo "Creating output and temporary directories" >> DLX_Ctrl_scenicplus_pipeline.log
mkdir -p outs
mkdir -p tmp

# Define paths to manually added files
GENOME_ANNOTATION_FILE="/faststorage/project/Hypothalamus/SCENIC/outs/genome_annotation_with_chr.tsv"
CHROMSIZES_FILE="/faststorage/project/Hypothalamus/SCENIC/outs/chromsizes.tsv"

echo "Using manually specified genome annotation file: $GENOME_ANNOTATION_FILE" >> DLX_Ctrl_scenicplus_pipeline.log
echo "Using manually specified chromsizes file: $CHROMSIZES_FILE" >> DLX_Ctrl_scenicplus_pipeline.log


# Modify the config.yaml file
echo "Modifying config.yaml file" >> DLX_Ctrl_scenicplus_pipeline.log
cat <<EOL > scplus_pipeline/Snakemake/config/config.yaml
input_data:
  cisTopic_obj_fname: "/faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT/DLX/CTRL_DLX_outs/CTRL_DLX_obj2.pkl"
  GEX_anndata_fname: "/faststorage/project/Hypothalamus/SCENIC/RNA/OUTPUT/DLX/CTRL_DLX2.h5ad"
  region_set_folder: "/faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT/DLX/CTRL_DLX_outs/region_sets"
  ctx_db_fname: "/faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT/DLX/CTRL_DLX_outs/CTRL_DLX_CTRL.regions_vs_motifs.rankings.feather"
  dem_db_fname: "/faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT/DLX/CTRL_DLX_outs/CTRL_DLX_CTRL.regions_vs_motifs.scores.feather"
  path_to_motif_annotations: "/faststorage/project/Hypothalamus/SCENIC/v10nr_clust_public/snapshots/motifs-v10-nr.mgi-m0.00001-o0.0.tbl"

output_data:
  combined_GEX_ACC_mudata: "/faststorage/project/Hypothalamus/SCENIC/outs/ACC_GEX.h5mu"
  dem_result_fname: "/faststorage/project/Hypothalamus/SCENIC/outs/dem_results.hdf5"
  ctx_result_fname: "/faststorage/project/Hypothalamus/SCENIC/outs/ctx_results.hdf5"
  output_fname_dem_html: "/faststorage/project/Hypothalamus/SCENIC/outs/dem_results.html"
  output_fname_ctx_html: "/faststorage/project/Hypothalamus/SCENIC/outs/ctx_results.html"
  cistromes_direct: "/faststorage/project/Hypothalamus/SCENIC/outs/cistromes_direct.h5ad"
  cistromes_extended: "/faststorage/project/Hypothalamus/SCENIC/outs/cistromes_extended.h5ad"
  tf_names: "/faststorage/project/Hypothalamus/SCENIC/outs/tf_names.txt"
  genome_annotation: "$GENOME_ANNOTATION_FILE"
  chromsizes: "$CHROMSIZES_FILE"
  search_space: "/faststorage/project/Hypothalamus/SCENIC/outs/search_space.tsv"
  tf_to_gene_adjacencies: "/faststorage/project/Hypothalamus/SCENIC/outs/tf_to_gene_adj.tsv"
  region_to_gene_adjacencies: "/faststorage/project/Hypothalamus/SCENIC/outs/region_to_gene_adj.tsv"
  eRegulons_direct: "/faststorage/project/Hypothalamus/SCENIC/outs/eRegulon_direct.tsv"
  eRegulons_extended: "/faststorage/project/Hypothalamus/SCENIC/outs/eRegulons_extended.tsv"
  AUCell_direct: "/faststorage/project/Hypothalamus/SCENIC/outs/AUCell_direct.h5mu"
  AUCell_extended: "/faststorage/project/Hypothalamus/SCENIC/outs/AUCell_extended.h5mu"
  scplus_mdata: "/faststorage/project/Hypothalamus/SCENIC/outs/scplusmdata.h5mu"

params_general:
  temp_dir: "/faststorage/project/Hypothalamus/SCENIC/tmp"
  n_cpu: 20
  seed: 666

params_data_preparation:
  bc_transform_func: "\"lambda x: x\""
  is_multiome: False  # Set to False for non-multiome data
  key_to_group_by: "Cluster_Pass2"
  nr_cells_per_metacells: 5
  direct_annotation: "Direct_annot"
  extended_annotation: "Orthology_annot"
  species: "mmusculus"  # Mouse species
  biomart_host: "nov2020.archive.ensembl.org"
  search_space_upstream: "1000 150000"
  search_space_downstream: "1000 150000"
  search_space_extend_tss: "10 10"

params_motif_enrichment:
  species: "mus_musculus"  # Mouse species
  annotation_version: "v10nr_clust"
  motif_similarity_fdr: 0.001
  orthologous_identity_threshold: 0.0
  annotations_to_use: "Direct_annot Orthology_annot"
  fraction_overlap_w_dem_database: 0.4
  dem_max_bg_regions: 500
  dem_balance_number_of_promoters: True
  dem_promoter_space: 1000
  dem_adj_pval_thr: 0.05
  dem_log2fc_thr: 1.0
  dem_mean_fg_thr: 0.0
  dem_motif_hit_thr: 3.0
  fraction_overlap_w_ctx_database: 0.4
  ctx_auc_threshold: 0.005
  ctx_nes_threshold: 3.0
  ctx_rank_threshold: 0.05

params_inference:
  tf_to_gene_importance_method: "GBM"
  region_to_gene_importance_method: "GBM"
  region_to_gene_correlation_method: "SR"
  order_regions_to_genes_by: "importance"
  order_TFs_to_genes_by: "importance"
  gsea_n_perm: 1000
  quantile_thresholds_region_to_gene: "0.85 0.90 0.95"
  top_n_regionTogenes_per_gene: "5 10 15"
  top_n_regionTogenes_per_region: ""
  min_regions_per_gene: 0
  rho_threshold: 0.05
  min_target_genes: 10
EOL

# Run the Snakemake pipeline
echo "Running Snakemake pipeline commands directly using CLI" >> DLX_Ctrl_scenicplus_pipeline.log
cd scplus_pipeline/Snakemake
snakemake --cores 20 --latency-wait 300

# Check if Snakemake pipeline ran successfully
if [ $? -ne 0 ]; then
    echo "Snakemake pipeline failed" >> DLX_Ctrl_scenicplus_pipeline.log
    exit 1
fi

# Print the end time of the job
echo "Job finished at: $(date)" >> DLX_Ctrl_scenicplus_pipeline.log