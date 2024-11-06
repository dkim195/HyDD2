#!/bin/bash
#SBATCH --account=Hypothalamus
#SBATCH -c 10
#SBATCH --mem=256g
#SBATCH --partition=normal
#SBATCH --time=8:00:00
#SBATCH --mail-type=END
#SBATCH --job-name=deeplab2
#SBATCH --mail-user=tkim@dandrite.au.dk

# Stop the script if any command fails
set -e

# Load Anaconda module and activate the Conda environment for Python
source /home/tkim/miniforge3/etc/profile.d/conda.sh
echo "Activating Conda environment: deeplabcut"
conda activate DEEPLABCUT

echo "Training..."
python /faststorage/project/Hypothalamus/script/deep/deeplabcut_training.py

echo "Evaluating..."
python /faststorage/project/Hypothalamus/script/deep/deeplabcut_evaluating.py

echo "Analyzing..."
python /faststorage/project/Hypothalamus/script/deep/deeplabcut_analyzing.py


