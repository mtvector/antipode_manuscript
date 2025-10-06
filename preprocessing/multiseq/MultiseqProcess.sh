#!/bin/bash
#$ -o ~/log
#$ -e ~/log
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l mem_free=96G
#$ -l h_rt=23:00:00

export PATH="~/utils/miniconda3/bin/:$PATH"
source ~/.bashrc
conda activate multiseq
function join_by { local IFS="$1"; shift; echo "$*"; }

MSOUT=/wynton/group/ye/mtschmitz/macaquedevbrain/MULTISEQmacaque/pollena-MULTIseq_E80
KOUT=/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto/E80-2019_Multi-seq_kOut
python /wynton/home/ye/mschmitz1/code/macaque-dev-brain/multiseq/MultiseqFuzzyWuzzy.py $MSOUT/MultiseqKey.txt $MSOUT/MULTIseq_E80_S97_L002_R1_001.fastq $MSOUT/../MultiseqIndices.txt $KOUT/cellranger_barcodes.tsv $KOUT 0 
IFS=$'\r\n' GLOBIGNORE='*' command eval  "XYZ=($(awk -F '\t' '{print $2}' $KOUT/multiseq_outs/features.tsv))"
VAR_LIST=$(join_by , "${XYZ[@]}")
rm $KOUT/multiseq_outs/*.gz
gzip $KOUT/multiseq_outs/*
rm -r $KOUT/multiseq_outs/GMMdemux/
#GMM-demux $KOUT/multiseq_outs/ "$VAR_LIST" -o $KOUT/multiseq_outs/GMMdemux/ --simplified $KOUT/multiseq_outs/GMMdemux/
python ~/code/macaque-dev-brain/multiseq/RunDemuxEM.py $KOUT


MSOUT=/wynton/group/ye/mtschmitz/macaquedevbrain/MULTISEQmacaque/pollena-MULTIseq_E90
KOUT=/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto/E90-2019_Multi-seq_kOut
python /wynton/home/ye/mschmitz1/code/macaque-dev-brain/multiseq/MultiseqFuzzyWuzzy.py $MSOUT/MultiseqKey.txt $MSOUT/MULTIseq_E90_S98_L002_R1_001.fastq $MSOUT/../MultiseqIndices.txt $KOUT/cellranger_barcodes.tsv $KOUT 0
IFS=$'\r\n' GLOBIGNORE='*' command eval  "XYZ=($(awk -F '\t' '{print $2}' $KOUT/multiseq_outs/features.tsv))"
VAR_LIST=$(join_by , "${XYZ[@]}")
rm $KOUT/multiseq_outs/*.gz
gzip $KOUT/multiseq_outs/*
rm -r $KOUT/multiseq_outs/GMMdemux/
#GMM-demux $KOUT/multiseq_outs/ "$VAR_LIST" -o $KOUT/multiseq_outs/GMMdemux/ --simplified $KOUT/multiseq_outs/GMMdemux/
python ~/code/macaque-dev-brain/multiseq/RunDemuxEM.py $KOUT


MSOUT=/wynton/group/ye/mtschmitz/macaquedevbrain/MULTISEQmacaque/pollena-MULTIseq_E65-1
KOUT=/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto/E65-2019A_AND_E65-2019B_Multi-seq_1_kOut
python /wynton/home/ye/mschmitz1/code/macaque-dev-brain/multiseq/MultiseqFuzzyWuzzy.py $MSOUT/MultiseqKey.txt $MSOUT/MULTIseq_E65-1_S99_L002_R1_001.fastq $MSOUT/../MultiseqIndices.txt $KOUT/cellranger_barcodes.tsv $KOUT 0
IFS=$'\r\n' GLOBIGNORE='*' command eval  "XYZ=($(awk -F '\t' '{print $2}' $KOUT/multiseq_outs/features.tsv))"
VAR_LIST=$(join_by , "${XYZ[@]}")
rm $KOUT/multiseq_outs/*.gz
gzip $KOUT/multiseq_outs/*
rm -r $KOUT/multiseq_outs/GMMdemux/
#GMM-demux $KOUT/multiseq_outs/ "$VAR_LIST" -o $KOUT/multiseq_outs/GMMdemux/ --simplified $KOUT/multiseq_outs/GMMdemux/
python ~/code/macaque-dev-brain/multiseq/RunDemuxEM.py $KOUT

MSOUT=/wynton/group/ye/mtschmitz/macaquedevbrain/MULTISEQmacaque/pollena-MULTIseq_E65-2
KOUT=/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto/E65-2019A_AND_E65-2019B_Multi-seq_2_kOut
python /wynton/home/ye/mschmitz1/code/macaque-dev-brain/multiseq/MultiseqFuzzyWuzzy.py $MSOUT/MultiseqKey.txt $MSOUT/MULTIseq_E65-2_S100_L002_R1_001.fastq $MSOUT/../MultiseqIndices.txt $KOUT/cellranger_barcodes.tsv $KOUT 0
IFS=$'\r\n' GLOBIGNORE='*' command eval  "XYZ=($(awk -F '\t' '{print $2}' $KOUT/multiseq_outs/features.tsv))"
VAR_LIST=$(join_by , "${XYZ[@]}")
rm $KOUT/multiseq_outs/*.gz
gzip $KOUT/multiseq_outs/*
rm -r $KOUT/multiseq_outs/GMMdemux/
#GMM-demux $KOUT/multiseq_outs/ "$VAR_LIST" -o $KOUT/multiseq_outs/GMMdemux/ --simplified $KOUT/multiseq_outs/GMMdemux/
python ~/code/macaque-dev-brain/multiseq/RunDemuxEM.py $KOUT

MSOUT=/wynton/group/ye/mtschmitz/macaquedevbrain/MULTISEQmacaque/pollena-MULTIseq_E65-3
KOUT=/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto/E65-2019A_AND_E65-2019B_Multi-seq_3_kOut
python /wynton/home/ye/mschmitz1/code/macaque-dev-brain/multiseq/MultiseqFuzzyWuzzy.py $MSOUT/MultiseqKey.txt $MSOUT/MULTIseq_E65-3_S101_L002_R1_001.fastq $MSOUT/../MultiseqIndices.txt $KOUT/cellranger_barcodes.tsv $KOUT 0
IFS=$'\r\n' GLOBIGNORE='*' command eval  "XYZ=($(awk -F '\t' '{print $2}' $KOUT/multiseq_outs/features.tsv))"
VAR_LIST=$(join_by , "${XYZ[@]}")
rm $KOUT/multiseq_outs/*.gz
gzip $KOUT/multiseq_outs/*
rm -r $KOUT/multiseq_outs/GMMdemux/
#GMM-demux $KOUT/multiseq_outs/ "$VAR_LIST" -o $KOUT/multiseq_outs/GMMdemux/ --simplified $KOUT/multiseq_outs/GMMdemux/
python ~/code/macaque-dev-brain/multiseq/RunDemuxEM.py $KOUT

