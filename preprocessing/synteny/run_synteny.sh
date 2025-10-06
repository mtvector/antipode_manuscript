#!/bin/bash
#SBATCH --job-name synteny
#SBATCH --output /home/matthew.schmitz/log/syn_%A.out # %A is replaced by job ID, %a by array index
#SBATCH --error /home/matthew.schmitz/log/syn_%A.err
#SBATCH --time 96:00:00
#SBATCH --partition celltypes
#SBATCH --mem 128gb
#SBATCH --ntasks 1

source ~/.bashrc

conda activate sequences

BLOCK_SIZE=10000
ANCHOR_SIZE=10000

mkdir -p /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/gene_lists/synteny/

/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/cactus-bin-v2.8.0/bin/halSynteny ~/Matthew/genome/hal/3_species/3g.hal --minBlockSize $BLOCK_SIZE --maxAnchorDistance $ANCHOR_SIZE /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/gene_lists/synteny/qQtM.psl --queryGenome Macaca_mulatta --targetGenome Mus_musculus
/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/cactus-bin-v2.8.0/bin/halSynteny ~/Matthew/genome/hal/3_species/3g.hal --minBlockSize $BLOCK_SIZE --maxAnchorDistance $ANCHOR_SIZE /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/gene_lists/synteny/qMtQ.psl --queryGenome Mus_musculus --targetGenome Macaca_mulatta
/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/cactus-bin-v2.8.0/bin/halSynteny ~/Matthew/genome/hal/3_species/3g.hal --minBlockSize $BLOCK_SIZE --maxAnchorDistance $ANCHOR_SIZE /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/gene_lists/synteny/qMtH.psl --queryGenome Mus_musculus --targetGenome Homo_sapiens
/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/cactus-bin-v2.8.0/bin/halSynteny ~/Matthew/genome/hal/3_species/3g.hal --minBlockSize $BLOCK_SIZE --maxAnchorDistance $ANCHOR_SIZE /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/gene_lists/synteny/qHtM.psl --queryGenome Homo_sapiens --targetGenome Mus_musculus
/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/cactus-bin-v2.8.0/bin/halSynteny ~/Matthew/genome/hal/3_species/3g.hal --minBlockSize $BLOCK_SIZE --maxAnchorDistance $ANCHOR_SIZE /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/gene_lists/synteny/qQtH.psl --queryGenome Macaca_mulatta --targetGenome Homo_sapiens
/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/cactus-bin-v2.8.0/bin/halSynteny ~/Matthew/genome/hal/3_species/3g.hal --minBlockSize $BLOCK_SIZE --maxAnchorDistance $ANCHOR_SIZE /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/gene_lists/synteny/qHtQ.psl --queryGenome Homo_sapiens --targetGenome Macaca_mulatta
