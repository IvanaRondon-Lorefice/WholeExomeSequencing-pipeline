#!/bin/bash
#SBATCH --job-name=03_SNVs_step1_PON
#SBATCH -o logs/03_SNVs_step1_PON.out
#SBATCH -e logs/03_SNVs_step1_PON.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=50:00:00
#SBATCH --partition=NORMAL
#SBATCH --mail-user=irondon@cicbiogune.es
#SBATCH --mail-type=FAIL

#source /vols/GPArkaitz_bigdata/irondon/AC-82_WESmouse_pipeline/00_conf_env.env
source /opt/ohpc/pub/apps/anaconda3/cic-env
conda activate /vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/conda_envs/WESenv

PON_DIR="/vols/GPArkaitz_bigdata/irondon/AC-82_WESmouse_pipeline/03_SNVs_GATK/03_SNVs_step1_PON/"
# Directory containing VCF files
#VCF_FILES=$(ls ${PON_DIRECTORY}*.vcf.gz | sed 's/^/-vcfs /' | tr '\n' ' ')

# Combine the normal calls using CreateSomaticPanelOfNormals.
gatk CreateSomaticPanelOfNormals \
   -vcfs ${PON_DIR}/Sample_26.vcf.gz \
   -vcfs ${PON_DIR}/Sample_27.vcf.gz \
   -vcfs ${PON_DIR}/Sample_28.vcf.gz \
   -vcfs ${PON_DIR}/Sample_29.vcf.gz \
   -vcfs ${PON_DIR}/Sample_30.vcf.gz \
   -vcfs ${PON_DIR}/Sample_33.vcf.gz \
   --min-sample-count 1 \
   -O ${PON_DIR}/PON/PON_min-sample-count-1_2.vcf.gz




