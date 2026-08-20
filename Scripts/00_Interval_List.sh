#!/bin/bash

#---------------------------------------------------------------------------------#
# Convert BED to Interval List Format for GATK's CollectReadCounts to quantify CNVs
# in WES data for mouse 

# Bed file was downloaded from: https://www.twistbioscience.com/resources/data-files/twist-mouse-exome-panel-bed-file
# for Twist Mouse Exome Panel

# For this library, we need to use the GRCm38. 
#---------------------------------------------------------------------------------#

source /opt/ohpc/pub/apps/anaconda3/cic-env
conda activate /vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/conda_envs/WESenv

picard BedToIntervalList I=/vols/GPArkaitz_bigdata/irondon/AC-82_WESmouse/Twist_Mouse_Exome_Target_Rev1_7APR20_Chr-Modified.bed \
    O=/vols/GPArkaitz_bigdata/irondon/AC-82_WESmouse/Twist_Mouse_Exome.interval_list \
    SD=/vols/GPArkaitz_bigdata/DATA_shared/Genomes_Rocky/Sequences/Mus_musculus.GRCm38.dna.primary_assembly.dict
