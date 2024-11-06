import os
import pandas as pd
import pickle
import numpy as np
import logging
from pycisTopic.cistopic_class import create_cistopic_object
from pycisTopic.lda_models import run_cgs_models_mallet, evaluate_models
from pycisTopic.clust_vis import run_umap
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features)
from pycisTopic.utils import region_names_to_coordinates

# Set up logging
logging.basicConfig(filename="pycistopic_neurogenesis.log", level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

logging.info("Script started.")

# Change to the appropriate directory
os.chdir("/faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT/neurogenesis")

logging.info(f"Changed to directory: {os.getcwd()}")

# Path to the fragment matrix
matrix_path = "neurogenesis_count_matrix.feather"

# Load the fragment matrix
if not os.path.exists(matrix_path):
    logging.error(f"The file '{matrix_path}' does not exist.")
    raise FileNotFoundError(f"The file '{matrix_path}' does not exist.")
fragment_matrix = pd.read_feather(matrix_path)
fragment_matrix.set_index('rownames', inplace=True)

# Create the cisTopicObject
cisTopic_obj = create_cistopic_object(fragment_matrix)

# Load cell annotations
cell_annotation_path = "neurogenesis_cell_annotation.tsv"
if not os.path.exists(cell_annotation_path):
    logging.error(f"The file '{cell_annotation_path}' does not exist.")
    raise FileNotFoundError(f"The file '{cell_annotation_path}' does not exist.")
cell_data = pd.read_table(cell_annotation_path, index_col=0)

# Add cell annotations to the cisTopicObject
cisTopic_obj.add_cell_data(cell_data)

out_dir = "neurogenesis_outs"
os.makedirs(out_dir, exist_ok=True)
pickle.dump(cisTopic_obj, open(os.path.join(out_dir, "neurogenesis_obj.pkl"), "wb"))

# Configure memory for MALLET
os.environ['MALLET_MEMORY'] = '600G'  # Increase MALLET memory allocation
mallet_path = "/faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT/Mallet-202108/bin/mallet"

# Run models
models = run_cgs_models_mallet(
    cisTopic_obj,
    n_topics=[50],
    n_cpu=30,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    mallet_path=mallet_path
)

pickle.dump(models, open(os.path.join(out_dir, "neurogenesis_models.pkl"), "wb"))

# Evaluate models
model = evaluate_models(
    models,
    select_model=50,
    return_model=True
)
cisTopic_obj.add_LDA_model(model)

# Run UMAP
run_umap(cisTopic_obj, target='cell', scale=True)

X_umap = pd.read_csv("neurogenesis_cell_embeddings.csv", index_col=0)
X_umap.columns = ["UMAP_1", "UMAP_2"]
X_umap.index = [idx + '___pycistopic' for idx in X_umap.index]
cisTopic_obj.projections['cell']['UMAP'] = X_umap

# Binarize topics
region_bin_topics_top_3k = binarize_topics(cisTopic_obj, method='ntop', ntop=3000, plot=True, num_columns=5)
region_bin_topics_otsu = binarize_topics(cisTopic_obj, method='otsu', plot=True, num_columns=5)
binarized_cell_topic = binarize_topics(cisTopic_obj, target='cell', method='li', plot=True, num_columns=5, nbins=100)

# Find differentially accessible regions (DARs)
imputed_acc_obj = impute_accessibility(cisTopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp=0.05,
    min_mean=0.0125,
    max_mean=3,
    max_disp=np.inf,
    n_bins=30,
    n_top_features=None,
    plot=False
)
# Log the length of variable_regions
logging.info(f'Number of variable regions: {len(variable_regions)}')

markers_dict = find_diff_features(
    cisTopic_obj,
    imputed_acc_obj,
    variable='Cluster_Pass2',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1),
    n_cpu=10,
    split_pattern='-'
)

logging.info("Number of DARs found:")
for x in markers_dict:
    logging.info(f"  {x}: {len(markers_dict[x])}")

pickle.dump(cisTopic_obj, open(os.path.join(out_dir, "neurogenesis_obj.pkl"), "wb"))

# Save results
os.makedirs(os.path.join(out_dir, "region_sets"), exist_ok=True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok=True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok=True)
os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok=True)

for topic in region_bin_topics_otsu:
    region_names_to_coordinates(region_bin_topics_otsu[topic].index).sort_values(["Chromosome", "Start", "End"]).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep="\t",
        header=False, index=False
    )
for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(region_bin_topics_top_3k[topic].index).sort_values(["Chromosome", "Start", "End"]).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep="\t",
        header=False, index=False
    )
for cell_type in markers_dict:
    region_names_to_coordinates(markers_dict[cell_type].index).sort_values(["Chromosome", "Start", "End"]).to_csv(
        os.path.join(out_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep="\t",
        header=False, index=False
    )

logging.info("Script finished.")