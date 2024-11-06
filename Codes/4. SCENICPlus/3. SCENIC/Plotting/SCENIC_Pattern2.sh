#!/bin/bash
#SBATCH --account=Hypothalamus
#SBATCH -c 4
#SBATCH --mem=128g
#SBATCH --partition=normal
#SBATCH --time=00:20:00
#SBATCH --mail-type=END
#SBATCH --job-name=pycistopic_PATTERN
#SBATCH --mail-user=tkim@dandrite.au.dk


# Stop the script if any command fails
set -e

# Source the Conda setup script to ensure conda commands are available in this script
echo "Sourcing Conda setup script..."
source /home/tkim/miniforge3/etc/profile.d/conda.sh

# Activate the Conda environment for Python
echo "Activating Conda environment: scenicplus"
conda activate scenicplus

# Run the Python scripts and capture the output and errors in log files
echo "Running Python script: PATTERN_PYCISTOPIC.py"
python /faststorage/project/Hypothalamus/script/scenicscript/SCENIC_Pattern2.py >> /faststorage/project/Hypothalamus/script/scenicscript/SCENIC_Pattern_Output.log 2>&1
echo "Python PATTERN script execution completed."

# Run the Python scripts and capture the output and errors in log files
echo "Running Python script: Neurogenesis_PYCISTOPIC.py"
python /faststorage/project/Hypothalamus/script/scenicscript/SCENIC_Neurogenesis2.py >> /faststorage/project/Hypothalamus/script/scenicscript/SCENIC_Neurogenesis_Output.log 2>&1
echo "Python Neurogenesis script execution completed."


echo "Script execution completed. SCENICPLUS preparation is completed"
