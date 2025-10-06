#!/bin/bash
#SBATCH --job-name solo
#SBATCH --output /home/matthew.schmitz/log/solo_%A.out # %A is replaced by job ID, %a by array index
#SBATCH --error /home/matthew.schmitz/log/solo_%A.err
#SBATCH --time 24:00:00
#SBATCH --partition celltypes
#SBATCH --gres=gpu:1
#SBATCH --mem 128gb
#SBATCH --ntasks 1

source ~/.bashrc

conda activate solo

python ~/Matthew/code/triple-dev-wb/preprocessing/DoubletSolo.py