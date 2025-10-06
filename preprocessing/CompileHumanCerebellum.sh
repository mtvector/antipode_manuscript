#!/bin/bash
#SBATCH --job-name cb
#SBATCH --output /home/matthew.schmitz/log/cb.out # %A is replaced by job ID, %a by array index
#SBATCH --error /home/matthew.schmitz/log/cb.err
#SBATCH --time 12:00:00
#SBATCH --partition celltypes
#SBATCH --mem 128gb
#SBATCH --ntasks 1

source ~/.bashrc
conda activate scanpy

# python ~/Matthew/code/triple-dev-wb/mdb/whole-brain/human/HumanKDCb1_Presupervision-cerebellum.py
python ~/Matthew/code/triple-dev-wb/mdb/whole-brain/mouse/MouseKDCb2_MouseCortexFillin.py