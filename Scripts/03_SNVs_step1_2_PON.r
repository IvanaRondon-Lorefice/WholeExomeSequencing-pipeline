#/opt/R/R-4.1.2/bin/R


################################################################################
###############################  DATA LIST  ####################################
################################################################################

# Configuration env file
conf_env <- "../00_conf_env.env"

# Loading necessary libraries and data
.libPaths("../R_libraries")
library("dotenv")
load_dot_env(file = conf_env)  # Carga las variables del archivo .env

# Setting based list directories
Base_dir <- Sys.getenv("BASE_DIR")
Data_shared <- Sys.getenv("DATA_SHARED")

# Only normal samples
samples <- strsplit(Sys.getenv("NORMAL_PON_SAMPLES"), ",")[[1]]

# Conda Activate
conda_activate <- Sys.getenv("CONDA_ACTIVATE")

# Conda environment WES
conda_env_wes <- paste0(Data_shared, Sys.getenv("CONDA_ENV_WES"))

# Mouse genome
genome <- paste0(Data_shared, Sys.getenv("GENOME"))

# Germline source
snps_indels_data <- paste0(Data_shared,Sys.getenv("SNPS_INDELS_DATA"))

# Mapped sequences
mapped_directory <- paste0(Base_dir, Sys.getenv("MAPPED_QC_DIRECTORY"))

# Output PON
pon_directory <- paste0(Base_dir,Sys.getenv("PON_DIRECTORY"))

# Interval list 
interval_list <- paste0(Base_dir, Sys.getenv("INTERVAL_LIST"))

# Jobs dir
jobs_dir <- paste0(Base_dir, Sys.getenv("JOBS_DIR"))

# Minimal samples counts
minim_sample_counts <-  Sys.getenv("MIN_SAMPLE_COUNTS_PON")

# Jobs requirements
cpu <- Sys.getenv("CPU_PON")
memory <- Sys.getenv("MEMORY_PON")
time <- Sys.getenv("TIME_PON")
partition <- Sys.getenv("PARTITION_PON")
mail <- Sys.getenv("MAIL")

# Name of the file (regarding the sample variable)
filename <- paste(jobs_dir,"03_SNVs_step1_2_PON.sh", sep = "") 

# File to save the variants detected
cat(
  "#!/bin/bash",
  "#SBATCH --job-name=03_SNVs_step1_2_PON",
  "#SBATCH -o logs/03_SNVs_step1_2_PON.out",
  "#SBATCH -e logs/03_SNVs_step1_2_PON.err",
  paste("#SBATCH --cpus-per-task=", cpu, sep = ""),
  paste("#SBATCH --mem=", memory, sep = ""),
  paste("#SBATCH --time=", time, sep = ""),
  paste("#SBATCH --partition=", partition, sep = ""),
  paste("#SBATCH --mail-user=",mail,sep = ""),
  "#SBATCH --mail-type=FAIL",

  "\n\n",
  paste("source ",conda_activate, sep = ""), 
  paste("conda activate ", conda_env_wes, sep = ""), 
  "\n",

  "VCF_FILES=\\'\\'",
  paste("for file in ", pon_directory, "*.vcf.gz; do", sep = ""),
  "  VCF_FILES=\\'$VCF_FILES -vcfs $file\\'",
  "done",
  "\n",

  paste("gatk CreateSomaticPanelOfNormals --min-sample-count 1 $VCF_FILES",
  " -O ", pon_directory,"PON/PON_min-sample-count-",minim_sample_counts,".vcf.gz", sep = ""),

  file = filename, sep = "\n", append = FALSE)

system(paste("sbatch", filename, sep = " "))

