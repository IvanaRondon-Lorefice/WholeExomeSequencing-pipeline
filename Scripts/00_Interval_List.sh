#!/bin/bash

#---------------------------------------------------------------------------------#
# Convert BED to Interval List Format for GATK's CollectReadCounts to quantify CNVs
# in WES data for mouse 

# Bed file was downloaded from: https://www.twistbioscience.com/resources/data-files/twist-mouse-exome-panel-bed-file
# for Twist Mouse Exome Panel

# For this library, we need to use the GRCm38. 
#---------------------------------------------------------------------------------#

source "$CONDA_ACTIVATE"
conda activate "$CONDA_ENV_WES"

picard BedToIntervalList \
    I="${BASE_DIR}/reference/Twist_Mouse_Exome_Target.bed" \
    O="${BASE_DIR}/reference/Twist_Mouse_Exome.interval_list" \
    SD="${DATA_SHARED}/${GENOME_DICT}"
