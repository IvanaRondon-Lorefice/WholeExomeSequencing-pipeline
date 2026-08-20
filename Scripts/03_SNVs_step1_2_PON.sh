#!/bin/bash
#SBATCH --job-name=03_SNVs_step1_PON
#SBATCH -o logs/03_SNVs_step1_PON.out
#SBATCH -e logs/03_SNVs_step1_PON.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=50:00:00
#SBATCH --partition=NORMAL
#SBATCH --mail-type=FAIL

source "$CONDA_ACTIVATE"
conda activate "$CONDA_ENV_WES"

PON_DIR="PON_dir/"
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
   --min-sample-count 1 \ # Optional
   -O ${PON_DIR}/PON_min-sample-count-1_2.vcf.gz




