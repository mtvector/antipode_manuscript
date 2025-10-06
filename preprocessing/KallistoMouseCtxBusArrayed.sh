#!/bin/bash
#SBATCH --export=ALL
#SBATCH --chdir=.
#SBATCH --error=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/extra_mouse_ctx/fastq/log/kalli-%A_%a.err
#SBATCH --output=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/extra_mouse_ctx/fastq/log/kalli-%A_%a.out
#SBATCH --job-name=kalli
#SBATCH --array=0-14
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=8:00:00
#SBATCH --partition=celltypes
# Optionally, request scratch space if your system supports it

DEST_DIR=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/extra_mouse_ctx/fastq
source ~/.bashrc
conda activate kallisto2

# Convert the Slurm array index (1-based) to a zero-based index.
TASK_ID=$((SLURM_ARRAY_TASK_ID - 1))
threads=8
genome=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/genome/kallisto/mm10

echo "Task ID: $TASK_ID"

# Build a list of unique sample prefixes by removing the _R1.fastq* suffix.
files=($(ls $DEST_DIR/*_R1.fastq* | sed -E 's/_R1\.fastq(.gz)?//' | sort | uniq | xargs -n 1 basename))
file="${files[$TASK_ID]}"
echo "Processing sample: $file"

cd $DEST_DIR
echo "Looking for output: ${file}_kOut/output.correct.sort.bus"

# Get the average read length from the first R1 file.
readlen=$(zcat ${file}_R1.fastq* | head -4 | awk '{if(NR%4==2){count++; bases += length}} END{print bases/count}')
echo "Read length: $readlen"


echo "Running: kallisto bus -t $threads -i ${genome}/cDNA_introns.idx -o ${file}_kOut/ -x 10xv2 ${file}_R*.fastq*"
kallisto bus -t $threads -i ${genome}/cDNA_introns.idx -o ${file}_kOut/ -x 10xv2 ${file}_R*.fastq*
cd ${file}_kOut/
mkdir -p cDNA_capture/ introns_capture/ spliced/ unspliced/ tmp/

# bustools correct -w /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/cellranger-8.0.1/lib/python/cellranger/barcodes/3M-february-2018.txt -p output.bus \
#   | bustools sort -T tmp/ -o output.correct.sort.bus -t $threads -
  
bustools whitelist -o ./WL.txt  output.bus
bustools correct -w ./WL.txt -p output.bus | bustools sort -T tmp/ -o output.correct.sort.bus -t $threads -


# mkdir -p ${file}_kOut/
# if [[ $readlen == 26 ]]; then
#     echo "Running: kallisto bus -t $threads -i ${genome}/cDNA_introns.idx -o ${file}_kOut/ -x 10xv2 ${file}_R*.fastq*"
#     kallisto bus -t $threads -i ${genome}/cDNA_introns.idx -o ${file}_kOut/ -x 10xv2 ${file}_R*.fastq*
#     cd ${file}_kOut/
#     mkdir -p cDNA_capture/ introns_capture/ spliced/ unspliced/ tmp/
#     bustools correct -w /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/cellranger-8.0.1/lib/python/cellranger/barcodes/3M-february-2018.txt -p output.bus \
#       | bustools sort -T tmp/ -o output.correct.sort.bus -t $threads -
# else
#     echo "Running: kallisto bus -t $threads -i ${genome}/cDNA_introns.idx -o ${file}_kOut/ -x 10xv3 ${file}_R*.fastq*"
#     kallisto bus -t $threads -i ${genome}/cDNA_introns.idx -o ${file}_kOut/ -x 10xv3 ${file}_R*.fastq*
#     cd ${file}_kOut/
#     mkdir -p cDNA_capture/ introns_capture/ spliced/ unspliced/ tmp/
#     bustools correct -w /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/cellranger-8.0.1/lib/python/cellranger/barcodes/3M-february-2018.txt -p output.bus \
#       | bustools sort -T tmp/ -o output.correct.sort.bus -t $threads -
# fi

cd ${file}_kOut/
mkdir -p all/ all_em/

bustools count -o all/a -g $genome/cDNA_introns_t2g.markintron.txt -e matrix.ec -t transcripts.txt --genecounts output.correct.sort.bus
bustools count -o all_em/a --em -g $genome/cDNA_introns_t2g.markintron.txt -e matrix.ec -t transcripts.txt --genecounts output.correct.sort.bus

# (Optional) If you wish to copy the output to another location, define JOB_DIR and uncomment the next line.
# cp -r $DEST_DIR/${file}_kOut $JOB_DIR
