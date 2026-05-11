#!/usr/bin/bash

# Run Cell Ranger multi for 10x Genomics Flex snucRNA-seq libraries from the postmenopausal human ovary atlas.
#
# Cell Ranger version used: v8.0.0
# Reference: GRCh38-2024-A
# Probe set: 10x Human Transcriptome Probe Set v1.0.1

# Define project directory
nostromo_datadir="/path/to/project_directory"
cd "${nostromo_datadir}/cellranger_results"

# Run Cell Ranger multi for each plex
for plex_num in {1..4}; do
  id="plex${plex_num}"
  configcsv="$nostromo_datadir/configs_multi_cellranger/plex${plex_num}_multiconfig.csv"

  cellranger multi --id $id --csv $configcsv
done

# Copy filtered per-sample .h5 count matrices to a central folder
h5_counts_dir="$nostromo_datadir/cellranger_results/count_matrices_h5"
for plex_num in {1..4}; do
    sample_outs_dir="$nostromo_datadir/cellranger_results/plex${plex_num}/outs/per_sample_outs"
    for subdir in $sample_outs_dir/*; do
        cp $subdir/count/sample_filtered_feature_bc_matrix.h5 $h5_counts_dir/$(basename $subdir)_P${plex_num}_filtered_counts.h5
    done
done
