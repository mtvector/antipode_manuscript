#!/bin/bash
#SBATCH --job-name=mapv1    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=256gb                     # Job memory request (per node)
#SBATCH --time=26:05:00               # Time limit hrs:min:sec
#SBATCH --output=map.log   # Standard output and error log
#SBATCH --partition celltypes         # Partition used for processing
#SBATCH --tmp=50G                     # Request the amount of space your jobs needs on /scratch/fast

source ~/.bashrc
conda activate mapcells

H5AD=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/mouse_wb/AIT21.all.freeze.230815_humanorthos.h5ad
REF_PATH=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/mouse_wb/ait_mapping_reference

python -m cell_type_mapper.cli.precompute_stats_scrattch \
--h5ad_path $H5AD \
--hierarchy '["class_label", "subclass_label", "supertype_label"]' \
--output_path $REF_PATH/precompute_stats.h5 \
--normalization log2CPM \
--tmp_dir ${REF_PATH}/temp/
echo "1 done"#,"cluster_label"

python -m cell_type_mapper.cli.reference_markers \
--precomputed_path_list '["'${REF_PATH}'/precompute_stats.h5"]' \
--output_dir ${REF_PATH}/ \
--tmp_dir ${REF_PATH}/temp/
echo "2 done"

QUERY_H5AD=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/extra/KDCbVelocityMouseWbAdultPresupervisionHumanOrthos.h5ad
OUT_PATH="${QUERY_H5AD%.h5ad}_ABC_MAPPING"
mkdir -p $OUT_PATH

python -m cell_type_mapper.cli.query_markers \
--reference_marker_path_list '["'${REF_PATH}'/reference_markers.h5"]' \
--output_path ${OUT_PATH}/query_markers.json
echo "4 done"

python -m cell_type_mapper.cli.from_specified_markers \
--query_path $QUERY_H5AD \
--query_markers.serialized_lookup ${OUT_PATH}/query_markers.json \
--type_assignment.normalization log2CPM \
--precomputed_stats.path $REF_PATH/precompute_stats.h5 \
--extended_result_path ${OUT_PATH}/hann_results.json \
> ${OUT_PATH}/log_outputs.txt 2>&1
echo "5 done"

python -m cell_type_mapper.cli.map_to_on_the_fly_markers \
--query_path $QUERY_H5AD \
--layer spliced \
--type_assignment.normalization raw \
--precomputed_stats.path $REF_PATH/precompute_stats.h5 \
--extended_result_path ${OUT_PATH}/hann_results.json \
> ${OUT_PATH}/log_outputs.txt 2>&1


