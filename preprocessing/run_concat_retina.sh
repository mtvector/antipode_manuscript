#!/bin/bash
#SBATCH --job-name concat_retina
#SBATCH --output /home/matthew.schmitz/log/concatr_%A.out # %A is replaced by job ID, %a by array index
#SBATCH --error /home/matthew.schmitz/log/concatr_%A.err
#SBATCH --time 24:00:00
#SBATCH --partition celltypes
#SBATCH --mem 511gb
#SBATCH --ntasks 1

source ~/.bashrc

conda activate scanpy

python ~/Matthew/code/triple-dev-wb/preprocessing/concat_retina.py