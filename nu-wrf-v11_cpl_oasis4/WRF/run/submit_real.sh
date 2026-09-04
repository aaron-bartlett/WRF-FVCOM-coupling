#!/bin/bash --login
#SBATCH -A glm200001 
#SBATCH --job-name=real_glm
#SBATCH --ntasks=16
#SBATCH --time=07:55:00
#SBATCH --output out_real_%j.txt
#SBATCH --error err_real_%j.txt

cd /compass/glm200001/cmu/coupled-run/nu-wrf-v11_cpl_oasis4/WRF/run

source /compass/glm200001/cmu/coupled-run/load_modules.sh

ln -s /compass/glm200001/cmu/coupled-run/nu-wrf-v11_cpl_oasis4/WPS/met_em* .

srun --verbose -n 16 --export=ALL /compass/glm200001/cmu/coupled-run/nu-wrf-v11_cpl_oasis4/WRF/run/real.exe 
