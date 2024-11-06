#!/bin/bash
#SBATCH --account=Hypothalamus
#SBATCH -c 30
#SBATCH --mem=700g
#SBATCH --partition=normal
#SBATCH --time=15:00:00
#SBATCH --mail-type=END
#SBATCH --job-name=pycistopic_neurogenesis
#SBATCH --mail-user=tkim@dandrite.au.dk


# Stop the script if any command fails
set -e

# Source the Conda setup script to ensure conda commands are available in this script
echo "Sourcing Conda setup script..."
source /home/tkim/miniforge3/etc/profile.d/conda.sh

# Activate the Conda environment for R
echo "Activating Conda environment: my_r_env"
conda activate my_r_env

# Run the R script and capture the output and errors in a log file
echo "Running R script: neurogenesis_SIGNAC2PYCISTOPIC.R"
Rscript /faststorage/project/Hypothalamus/script/scenicscript/NEUROGENESIS_SIGNAC2PYCISTOPIC.R > /faststorage/project/Hypothalamus/SCENIC/NEUROGENESIS_SIGNAC2PYCISTOPIC.log 2>&1
echo "R script execution completed."

# Deactivate Conda environment
echo "Deactivating Conda environment: my_r_env"
conda deactivate

# Activate the Conda environment for Python
echo "Activating Conda environment: scenicplus"
conda activate scenicplus

# Run the Python scripts and capture the output and errors in log files
echo "Running Python script: neurogenesis_PYCISTOPIC.py"
python /faststorage/project/Hypothalamus/script/scenicscript/NEUROGENESIS_PYCISTOPIC.py >> /faststorage/project/Hypothalamus/SCENIC/NEUROGENESIS_SIGNAC2PYCISTOPIC.log 2>&1
echo "Python neurogenesis script execution completed."


echo "Script execution completed. SCENICPLUS neurogenesis preparation is completed"
