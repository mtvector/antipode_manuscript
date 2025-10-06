#!/bin/bash
#SBATCH --export=ALL
#SBATCH --error=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/extra_mouse_ctx/fastq/log/cbk-%A_%a.err
#SBATCH --output=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/extra_mouse_ctx/fastq/log/cbk-%A_%a.out
#SBATCH --job-name=cellbender_remove_bg
#SBATCH --partition celltypes         # Partition used for processing
#SBATCH --array=1-14
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=9:59:59
#SBATCH --gres=gpu:1 --constraint="a100|v100"#

JOB_DIR=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/extra_mouse_ctx/fastq

source ~/.bashrc
conda activate cellbender

cd $JOB_DIR
# Slurm automatically sets CUDA_VISIBLE_DEVICES for allocated GPUs,
# but you can override it if needed:
# export CUDA_VISIBLE_DEVICES=$SLURM_JOB_GPUS

# Adjust array index (bash arrays are 0-indexed)
TASK_ID=$((SLURM_ARRAY_TASK_ID - 1))
EPOCH=200
MAXCELLS=8000
MINCELLS=$MAXCELLS
ADDCELLS=21000
FILES=($JOB_DIR/*_kOut)
FILE="${FILES[$TASK_ID]}"
FILEBASE=${FILE##*/}

# Determine the number of cells from the emptydrops file
CELLS=$(wc -l < "$FILE/all_em/barcodes_emptydrops.tsv")
if [ -z "$CELLS" ]; then 
    CELLS=$MAXCELLS
fi
MINCELLS=$(($CELLS < $MAXCELLS ? $CELLS : $MAXCELLS))
DROPS=$((MINCELLS + ADDCELLS))

# MINCOUNTS=$(($CELLS <= $MAXCELLS ? 20 : 100))
# MINCOUNTS=10

DIRNAME=aem_cellbended_150_750_200e_V0.2
echo "$FILEBASE"
if [ -r "$FILE/$DIRNAME/aem_cellbended.pdf" ]; then
    echo "Already processed."
else
    cd "$FILE"
    if [ -r "$FILE/$DIRNAME/aem_cellbended.log" ]; then
        echo "Log file found."
    else
        echo "No log file found."
    fi
    echo "Drops: $DROPS"
    echo "Cells: $CELLS"
    echo "Epochs: $EPOCH"
    mkdir -p ./"$DIRNAME"/
    cellbender remove-background \
                 --input ./all_em/ \
                 --output ./"$DIRNAME"/aem_cellbended.h5 \
                 --expected-cells "$MINCELLS" \
                 --total-droplets-included "$DROPS" \
                 --epochs "$EPOCH" \
                 --z-dim 150 --z-layers 750 \
                 --empty-drop-training-fraction .7 \
                 --model ambient \
                 --low-count-threshold 10 \
                 --cuda

                # --low-count-threshold "$MINCOUNTS" \

fi
exit 1
