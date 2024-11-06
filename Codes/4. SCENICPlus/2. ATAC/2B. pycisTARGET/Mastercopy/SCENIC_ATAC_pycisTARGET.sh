#!/bin/bash
#SBATCH --account=Hypothalamus
#SBATCH -c 20
#SBATCH --mem=512g
#SBATCH --partition=normal
#SBATCH --time=15:00:00
#SBATCH --mail-type=END
#SBATCH --job-name=CISTARGET
#SBATCH --mail-user=tkim@dandrite.au.dk


#OBJECT TO FILE NAME
# Stop the script if any command fails
set -e

# Source the Conda setup script to ensure conda commands are available in this script
echo "Sourcing Conda setup script..."
source /home/tkim/miniforge3/etc/profile.d/conda.sh

# Activate the Conda environment for Python
echo "Activating Conda environment: scenicplus"
conda activate scenicplus

# Run the Python scripts and capture the output and errors in log files
echo "Running Python script: DLX_CTRL_Consensus.py"
python /faststorage/project/Hypothalamus/script/scenicscript/OBJECT_Consensus.py >> /faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT/OBJECT_Consensus.log 2>&1

# Deactivate Conda environment
echo "Deactivating Conda environment: my_r_env"
conda deactivate

# Activate the Conda environment for Python
echo "Activating Conda environment: create_cisTarget_databases"
conda activate create_cisTarget_databases

# Set PATH explicitly to include the directory containing bedtools
export PATH="$HOME/miniforge3/envs/create_cisTarget_databases/bin:$PATH"

# Debugging: Check if bedtools is accessible
echo "Checking bedtools path:"
which bedtools

# Define common variables
GENOME_FASTA="/home/tkim/Hypothalamus/SCENIC/mm10/mm10.fa"
CHROMSIZES="/home/tkim/Hypothalamus/SCENIC/mm10/mm10.chrom.sizes"
DATABASE_PREFIX="OBJECT"
SCRIPT_DIR="/home/tkim/Hypothalamus/script/scenicscript/create_cisTarget_databases"
CBDIR="/home/tkim/Hypothalamus/SCENIC/v10nr_clust_public/singletons"
MOTIF_LIST="/home/tkim/Hypothalamus/SCENIC/v10nr_clust_public/singletons/motifs.txt"
REGION_BED="/faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT/DLX/OBJECT_outs/consensus_peak_calling/consensus_regions.bed"
CTRL_OUT_DIR="/faststorage/project/Hypothalamus/SCENIC/ATAC/OUTPUT/DLX/OBJECT_outs"
FASTA_FILE="${CTRL_OUT_DIR}/mm10.OBJECT.fa"

# Ensure output directory exists
mkdir -p ${CTRL_OUT_DIR}

echo "Processing CTRL DLx..."
echo "Creating FASTA for CTRL without padding..."
${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
    ${GENOME_FASTA} \
    ${CHROMSIZES} \
    ${REGION_BED} \
    ${FASTA_FILE} \
    1000 \
    yes >> ${CTRL_OUT_DIR}/create_fasta.log 2>&1

echo "Listing motifs..."
ls ${CBDIR} > ${MOTIF_LIST}

#had to manually make motif.txt.cb file from txt file using ls
#this step takes forever, 12 hrs for this set for ctrl_dlx

echo "Creating cistarget motif databases for CTRL..."
${SCRIPT_DIR}/create_cistarget_motif_databases.py \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${CTRL_OUT_DIR}/${DATABASE_PREFIX}_CTRL \
	--bgpadding 1000 \
    -t 20 >> ${CTRL_OUT_DIR}/create_cistarget_db.log 2>&1

echo "Script execution completed."