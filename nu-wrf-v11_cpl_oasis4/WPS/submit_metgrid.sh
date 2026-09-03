#!/bin/bash --login
#SBATCH -A glm200001
#SBATCH --job-name=metgrid_glm
#SBATCH --ntasks=8
#SBATCH --time=08:00:00

cd /compass/glm200001/cmu/coupled-run/nu-wrf-v11_cpl_oasis4/WPS

source /compass/glm200001/cmu/coupled-run/load_modules.sh

srun --verbose -n 8 --export=ALL /compass/glm200001/bart753/coupled_run/WPS/metgrid.exe

