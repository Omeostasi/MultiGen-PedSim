#!/bin/bash
#SBATCH --job-name=blupf90
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --time=72:00:00
#SBATCH --account=arch_sim

# Initialize conda for bash
source /home/giovincavallo/miniforge3/etc/profile.d/conda.sh

# Activate environment
conda activate blupf90

# Generate real values from renumf90 par file (the ouput of the previous step)

gibbsf90+ renf90.par --samples 10000 --burnin 5000 --interval 10  

# follows postgibbs
