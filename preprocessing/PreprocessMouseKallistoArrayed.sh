#!/bin/bash
#SBATCH --export=ALL
#SBATCH --error=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/extra_mouse_ctx/fastq/log/kalli-%A_%a.err
#SBATCH --output=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/extra_mouse_ctx/fastq/log/kalli-%A_%a.out
#SBATCH --job-name=kalli
#SBATCH --array=1-15
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --time=2:59:00
#SBATCH --partition=celltypes

JOB_DIR=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/extra_mouse_ctx/fastq
GENOME=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/genome/kallisto/mm10

source "/wynton/home/ye/mschmitz1/.bashrc"
conda activate kallisto2

cd $JOB_DIR
export CUDA_VISIBLE_DEVICES=$SLURM_JOB_GPUS
TASK_ID=$((SLURM_ARRAY_TASK_ID - 1))
FILES=($JOB_DIR/*_kOut)
FILE="${FILES[$TASK_ID]}"
FILEBASE=${FILE##*/}

cd $FILE/all_em
rm -f matrix.mtx.gz
rm -f features.tsv.gz
rm -f barcodes.tsv.gz
cp a.barcodes.txt barcodes.tsv

# Write Matrix Market header, then reformat a.mtx and append to matrix.mtx.
echo "%%MatrixMarket matrix coordinate real general\n%\n%" > matrix.mtx
sed 's/%.*//' a.mtx | awk -F" " '{print $2, $1, $3}' >> matrix.mtx

# Create features.tsv using information from the genome file and a.genes.txt.
awk -F" " 'NR==FNR { a[$2]=$2; b[$2]=$3; next } { print a[$1] "\t" b[$1] "\tGene Expression" }' \
    $GENOME/cDNA_introns_t2g.markintron.txt a.genes.txt > features.tsv

gzip *.tsv
gzip matrix.mtx
