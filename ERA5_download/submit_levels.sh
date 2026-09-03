#!/bin/bash --login
#SBATCH -A glm200001 
#SBATCH --job-name=dl_era5plv
#SBATCH --partition=normal
#SBATCH --time=08:00:00               
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=aaronbar@umich.edu
cd /compass/glm200001/bart753/ERA5-download
module load py-pip
source /home/bart753/.venv/bin/activate
python3 wget_levels.py > log.download_lev 2>&1

