#!/bin/bash
#SBATCH --job-name nbconvert_array
#SBATCH --output /home/matthew.schmitz/log/nbconvert_%A_%a.out # %A is replaced by job ID, %a by array index
#SBATCH --error /home/matthew.schmitz/log/nbconvert_%A_%a.err
#SBATCH --time 225:00:00
#SBATCH --partition celltypes
#SBATCH --array=135# Create an array job for each line in notebooks.txt
#SBATCH --gres=gpu:1 --constraint="a100"#|v100
#SBATCH --mem 148gb
#SBATCH --ntasks 1
#SBATCH --begin=now#+1hour

source ~/.bashrc

conda activate antipode

# Read the notebook path for the current array task
NOTEBOOK=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/Matthew/code/triple-dev-wb/model_runs/notebooks.txt)

jupyter nbconvert --ExecutePreprocessor.allow_errors=True --to notebook --execute "${NOTEBOOK}" --inplace
