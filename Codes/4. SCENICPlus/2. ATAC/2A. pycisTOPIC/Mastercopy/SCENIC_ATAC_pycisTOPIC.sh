#!/bin/bash
#SBATCH --account=Hypothalamus
#SBATCH -c 10
#SBATCH --mem=512g
#SBATCH --partition=normal
#SBATCH --time=8:00:00
#SBATCH --mail-type=END
#SBATCH --job-name=pycistopic_OBJECT
#SBATCH --mail-user=tkim@dandrite.au.dk

#Change OBJECT TO THE FILE NAME

# Stop the script if any command fails
set -e

# Source the Conda setup script to ensure conda commands are available in this script
echo "Sourcing Conda setup script..."
source /home/tkim/miniforge3/etc/profile.d/conda.sh

# Activate the Conda environment for R
echo "Activating Conda environment: my_r_env"
conda activate my_r_env

# Run the R script and capture the output and errors in a log file
echo "Running R script: OBJECT_SIGNAC2PYCISTOPIC.R"
Rscript /faststorage/project/Hypothalamus/script/scenicscript/OBJECT_SIGNAC2PYCISTOPIC.R > /faststorage/project/Hypothalamus/SCENIC/OBJECT_SIGNAC2PYCISTOPIC.log 2>&1
echo "R script execution completed."

# Deactivate Conda environment
echo "Deactivating Conda environment: my_r_env"
conda deactivate

# Activate the Conda environment for Python
echo "Activating Conda environment: scenicplus"
conda activate scenicplus

# Run the Python scripts and capture the output and errors in log files
echo "Running Python script: OBJECT_PYCISTOPIC.py"
python /faststorage/project/Hypothalamus/script/scenicscript/OBJECT_PYCISTOPIC.py >> /faststorage/project/Hypothalamus/SCENIC/OBJECT_SIGNAC2PYCISTOPIC.log 2>&1
echo "Python script execution completed."

echo "Script execution completed. SCENICPLUS OBJECT preparation is completed"
