#!/bin/bash
#$ -o ~/log
#$ -e ~/log
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l mem_free=96G
#$ -l h_rt=1:00:00
#$ -t 1-5:1

export PATH="~/utils/miniconda3/bin/:$PATH"
source ~/.bashrc
source ~/.bash_profile

source activate multiseq
conda activate multiseq


JOB_DIR=/wynton/scratch/mtschmitz/sara_data
cd $JOB_DIR

SGE_TASK_ID=`expr $SGE_TASK_ID - 1`
files=($(ls -d $JOB_DIR/*_mOut | awk -F'_mOut' '{print $1}' | uniq))
file="${files[$SGE_TASK_ID]}"
#file="${files[3]}"
echo "FILE"
echo $file

KOUT=${file}_mOut
cd $KOUT
fqfile=($(ls $KOUT/MULTISEQ*.fastq*))
fqfile="${fqfile[0]}"
echo $fqfile
python /wynton/home/ye/mschmitz1/code/macaque-dev-brain/multiseq/MultiseqFuzzyWuzzyGeneral.py $JOB_DIR/MultiseqSamples.txt $fqfile ~/code/macaque-dev-brain/multiseq/MultiseqIndices.txt $KOUT/outs/filtered_feature_bc_matrix/barcodes.tsv.gz $KOUT 0 

python ~/code/macaque-dev-brain/multiseq/RunDemuxEM.py $KOUT outs/filtered_feature_bc_matrix

function join_by { local IFS="$1"; shift; echo "$*"; }
IFS=$'\r\n' GLOBIGNORE='*' command eval  "XYZ=($(awk -F '\t' '{print $2}' $KOUT/multiseq_outs/features.tsv))"
VAR_LIST=$(join_by , "${XYZ[@]}")
rm $KOUT/multiseq_outs/*.gz
gzip $KOUT/multiseq_outs/*
rm -r $KOUT/multiseq_outs/GMMdemux/
GMM-demux $KOUT/multiseq_outs/ "$VAR_LIST" -o $KOUT/multiseq_outs/GMMdemux/ --simplified $KOUT/multiseq_outs/GMMdemux/
