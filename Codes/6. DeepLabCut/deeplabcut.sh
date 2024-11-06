#!/bin/bash
#SBATCH --account=Hypothalamus
#SBATCH -c 10
#SBATCH --mem=256g
#SBATCH --partition=normal
#SBATCH --time=12:00:00
#SBATCH --mail-type=END
#SBATCH --job-name=deeplab1
#SBATCH --mail-user=tkim@dandrite.au.dk

# Stop the script if any command fails
set -e

# Load Anaconda module and activate the Conda environment for Python
source /home/tkim/miniforge3/etc/profile.d/conda.sh
echo "Activating Conda environment: deeplabcut"
conda activate DEEPLABCUT

# Create Project and Extract Frames
echo "Creating Project and Extracting Frames..."
python /faststorage/project/Hypothalamus/script/deep/deeplabcut_create_extract.py

echo "Proceed to next step..."

#need to manually annotate/label frames
#need to check annoated frames
