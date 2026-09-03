
#!/bin/bash --login
#SBATCH -A glm200001
#SBATCH --job-name=coupled_run_glm
#SBATCH --ntasks=125                # must equal WRF_TASKS + FVCOM_TASKS below
#SBATCH --time=01:00:00             # short test window -- extend for production
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

...
