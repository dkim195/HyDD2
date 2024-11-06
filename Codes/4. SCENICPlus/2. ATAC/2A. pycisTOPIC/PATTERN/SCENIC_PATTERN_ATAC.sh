#!/bin/bash
#SBATCH --account=Hypothalamus
#SBATCH -c 20
#SBATCH --mem=512g
#SBATCH --partition=normal
#SBATCH --time=15:00:00
#SBATCH --mail-type=END
#SBATCH --job-name=pycistopic_PATTERN
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
echo "Running R script: PATTERN_SIGNAC2PYCISTOPIC.R"
Rscript /faststorage/project/Hypothalamus/script/scenicscript/PATTERN_SIGNAC2PYCISTOPIC.R > /faststorage/project/Hypothalamus/SCENIC/PATTERN_SIGNAC2PYCISTOPIC.log 2>&1
echo "R script execution completed."

# Deactivate Conda environment
echo "Deactivating Conda environment: my_r_env"
conda deactivate

# Activate the Conda environment for Python
echo "Activating Conda environment: scenicplus"
conda activate scenicplus

# Run the Python scripts and capture the output and errors in log files
echo "Running Python script: PATTERN_PYCISTOPIC.py"
python /faststorage/project/Hypothalamus/script/scenicscript/PATTERN_PYCISTOPIC.py >> /faststorage/project/Hypothalamus/SCENIC/PATTERN_SIGNAC2PYCISTOPIC.log 2>&1
echo "Python PATTERN script execution completed."


echo "Script execution completed. SCENICPLUS PATTERN preparation is completed"
