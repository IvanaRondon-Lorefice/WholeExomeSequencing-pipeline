#!/bin/bash

# Warning: Make sure the 00_conf_env.env is well formated. If needed run dos2unix

# Loading the configuration space
CONF_ENV="00_conf_env.env"

source ${CONF_ENV}

# Creating required folders
mkdir -p ${BASE_DIR}Scripts/jobs  
mkdir -p ${BASE_DIR}Scripts/logs

# Preprocessing
mkdir -p ${BASE_DIR}00_Trimming
mkdir -p ${BASE_DIR}01_FASTCs

# Mapping
mkdir -p ${BASE_DIR}02_Mapping
mkdir -p ${BASE_DIR}02_Mapping/02_Mapped_Sequences_QC
mkdir -p ${BASE_DIR}02_Mapping/Recal_Data_BQSR
mkdir -p ${BASE_DIR}02_Mapping/Marked_Duplicate_Metrics
mkdir -p ${BASE_DIR}02_Mapping/02_Mapped_Sequences

# SNVs analysis with GATK
mkdir -p ${BASE_DIR}/03_SNVs_GATK

# Conditionally create PON-related folders
if [ "${PON}" == "Yes" ]; then
    mkdir -p ${BASE_DIR}03_SNVs_GATK/03_SNVs_step1_PON
    mkdir -p ${BASE_DIR}03_SNVs_GATK/03_SNVs_step1_PON/PON
    mkdir -p ${BASE_DIR}03_SNVs_GATK/03_SNVs_step2_TUMOR-PON
    mkdir -p ${BASE_DIR}03_SNVs_GATK/03_SNVs_step2_TUMOR-PON/03_SNVs_step2_1_TUMOR_METRICS_BCFTOOLS
fi

# Conditionally create MATCH-related folders
if [ "${MATCH}" == "Yes" ]; then
    mkdir -p ${BASE_DIR}03_SNVs_GATK/3_SNVs_step2_TUMOR-MATCH
    mkdir -p ${BASE_DIR}03_SNVs_GATK/3_SNVs_step2_TUMOR-MATCH/03_SNVs_step2_1_TUMOR_METRICS_BCFTOOLS
fi
