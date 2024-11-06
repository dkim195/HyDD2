import os
import pandas as pd
import numpy as np
import logging
import subprocess
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Function to check if directory exists and create if not
def create_dir(directory):
    try:
        os.makedirs(directory, exist_ok=True)
        logging.info(f"Directory created or already exists: {directory}")
    except Exception as e:
        logging.error(f"Failed to create directory {directory}: {e}")
        raise

# Function to sanitize filenames to avoid issues with special characters
def sanitize_filename(filename):
    # Remove parentheses, replace spaces with underscores, and strip unwanted characters
    return filename.replace('(', '').replace(')', '').replace(' ', '_').replace(',', '')

# Set working directory and output directory
try:
    os.chdir("/faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT/neurogenesis")
    out_dir = "neurogenesis_outs"
    create_dir(out_dir)
except Exception as e:
    logging.error(f"Failed to set working directory or create output directory: {e}")
    raise

# Dictionary for fragment paths
fragments_dict = {
    "TK76": "/faststorage/project/Hypothalamus/data/scATAC/Fragments/TK76/fragments.tsv.gz",
    "TK77": "/faststorage/project/Hypothalamus/data/scATAC/Fragments/TK77/fragments.tsv.gz",
    "TK53": "/faststorage/project/Hypothalamus/data/scATAC/Fragments/TK53/fragments.tsv.gz",
    "TK54": "/faststorage/project/Hypothalamus/data/scATAC/Fragments/TK54/fragments.tsv.gz",
    "TK78": "/faststorage/project/Hypothalamus/data/scATAC/Fragments/TK78/fragments.tsv.gz",
    "TK79": "/faststorage/project/Hypothalamus/data/scATAC/Fragments/TK79/fragments.tsv.gz",
    "TK86": "/faststorage/project/Hypothalamus/data/scATAC/Fragments/TK86/fragments.tsv.gz",
    "TK87": "/faststorage/project/Hypothalamus/data/scATAC/Fragments/TK87/fragments.tsv.gz",
    "TK88": "/faststorage/project/Hypothalamus/data/scATAC/Fragments/TK88/fragments.tsv.gz",
    "TK89": "/faststorage/project/Hypothalamus/data/scATAC/Fragments/TK89/fragments.tsv.gz"

}
logging.info("Fragments dictionary defined")

# Check and load cell annotation data
cell_data_path = "/faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT/neurogenesis/neurogenesis_cell_annotation.tsv"
if not os.path.exists(cell_data_path):
    logging.error(f"Cell annotation data not found: {cell_data_path}")
    raise FileNotFoundError(f"Cell annotation data not found: {cell_data_path}")

try:
    cell_data = pd.read_table(cell_data_path, index_col=0)
    logging.info("Cell annotation data loaded")
except Exception as e:
    logging.error(f"Failed to load cell annotation data from {cell_data_path}: {e}")
    raise

# Check and load chromosome sizes
chromsizes_url = "https://hgdownload-test.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
try:
    chromsizes = pd.read_table(chromsizes_url, header=None, names=["Chromosome", "End"])
    chromsizes.insert(1, "Start", 0)
    logging.info("Chromosome sizes loaded")
except Exception as e:
    logging.error(f"Failed to load chromosome sizes from {chromsizes_url}: {e}")
    raise

# Create necessary directories for peak calling outputs
peak_calling_dirs = ["consensus_peak_calling", "consensus_peak_calling/pseudobulk_bed_files", "consensus_peak_calling/pseudobulk_bw_files"]
for directory in peak_calling_dirs:
    create_dir(os.path.join(out_dir, directory))

# Ensure correct barcode format and mapping to paths
try:
    cell_data.index = cell_data.index.map(lambda barcode: '_'.join(barcode.split('_')[1:] + [barcode.split('_')[0]]))
    cell_data['barcode_sample'] = cell_data.index.map(lambda x: x.split('_')[-1].replace('_', '').replace('-', ''))
    cell_data['barcode_sample'] = cell_data['barcode_sample'].astype(str)
    cell_data['fragments_path'] = cell_data['barcode_sample'].map(fragments_dict)
    fragments_path_dict = {sample_id: fragments_dict[sample_id] for sample_id in cell_data['barcode_sample'].unique()}
    logging.info("Barcode format transformed and fragments path mapped")
except Exception as e:
    logging.error(f"Failed to transform barcode format or map fragments paths: {e}")
    raise

# Sanitize the Clusters column specifically
cell_data['Cluster_Pass2'] = cell_data['Cluster_Pass2'].apply(sanitize_filename)

# Process pseudobulk data
try:
    bw_paths, bed_paths = export_pseudobulk(
        input_data=cell_data,
        variable="Cluster_Pass2",
        sample_id_col="barcode_sample",
        chromsizes=chromsizes,
        bed_path=os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
        bigwig_path=os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"),
        path_to_fragments=fragments_path_dict,
        n_cpu=20,
        normalize_bigwig=True,
        split_pattern="_"
    )
    logging.info("Pseudobulk peak calling completed")
except Exception as e:
    logging.error(f"Failed to execute pseudobulk peak calling: {e}")
    raise

# Save paths
def save_paths(file_path, paths_dict):
    try:
        with open(file_path, "w") as file:
            for key, value in paths_dict.items():
                file.write(f"{key}\t{value}\n")
        logging.info(f"Paths saved to {file_path}")
    except Exception as e:
        logging.error(f"Failed to save paths to {file_path}: {e}")
        raise

save_paths(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv"), bw_paths)
save_paths(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), bed_paths)

# Load paths
def load_paths(file_path):
    try:
        with open(file_path) as file:
            return dict(line.strip().split("\t") for line in file)
    except Exception as e:
        logging.error(f"Failed to load paths from {file_path}: {e}")
        raise

bw_paths = load_paths(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv"))
bed_paths = load_paths(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"))


# Peak calling
try:
    narrow_peak_dict = peak_calling(
        macs_path="macs2",
        bed_paths=bed_paths,
        outdir=os.path.join(out_dir, "consensus_peak_calling/MACS"),
        genome_size='mm',
        n_cpu=20,
        input_format='BEDPE',
        shift=73,
        ext_size=146,
        keep_dup='all',
        q_value=0.05
    )
    logging.info("MACS2 peak calling completed")
except Exception as e:
    logging.error(f"Failed during MACS2 peak calling: {e}")
    raise

# Get consensus peaks
try:
    consensus_peaks = get_consensus_peaks(
        narrow_peaks_dict=narrow_peak_dict,
        peak_half_width=250,
        chromsizes=chromsizes,
        path_to_blacklist="/home/tkim/Hypothalamus/SCENIC/mm10-blacklist.v2.bed"
    )
    logging.info("Consensus peaks generated")
except Exception as e:
    logging.error(f"Failed to generate consensus peaks: {e}")
    raise

# Save consensus peaks
try:
    consensus_peaks.to_bed(
        path=os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed"),
        keep=True,
        compression='infer',
        chain=False
    )
    logging.info("Consensus peaks saved to BED file")
except Exception as e:
    logging.error(f"Failed to save consensus peaks to BED file: {e}")
    raise