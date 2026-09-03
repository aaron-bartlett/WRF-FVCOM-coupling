#!/bin/bash --login
#SBATCH -A glm200001 
#SBATCH --job-name=dl_era5sfc
#SBATCH --partition=normal
#SBATCH --time=08:00:00               
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=aaronbar@umich.edu
cd /compass/glm200001/bart753/ERA5-download
module load py-pip
source /home/bart753/.venv/bin/activate
python3 cdsapi-surface.py > log.download_surf 2>&1
