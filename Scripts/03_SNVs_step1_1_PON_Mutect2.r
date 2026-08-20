#/opt/R/R-4.1.2/bin/R

################################################################################
####### THIS CODES CREATES A SCRIPT FOR EACH FILE FASTA AND SUMIT JOBS #########
################################################################################

#----------------------------------- MUTECT2 ----------------------------------

# Variant calling is the process of identifying genetic variations, such as single nucleotide polymorphisms (SNPs), 
# insertions, deletions, and structural variants, from DNA sequencing data. It is a crucial step in genomic analysis, 
# providing insights into genetic differences between individuals and populations. Mutect2, developed by the Broad Institute, 
# is a widely used tool for somatic variant calling in cancer genomics. It is specifically designed to identify somatic 
# mutations in tumor-normal paired samples with high sensitivity and specificity, making it ideal for detecting low-frequency 
# variants in complex tumor genomes. Mutect2's sophisticated algorithms and advanced filtering strategies enable accurate 
# identification of somatic mutations, helping researchers and clinicians unravel the genetic basis of cancer and guide 
# personalized treatment decisions.

#------------------------------------------------------------------------------

################################################################################


################################################################################
###############################  DATA LIST  ####################################
################################################################################

# Configuration env file
conf_env <- "../00_conf_env.env"

# Loading necessary libraries and data
.libPaths("/vols/GPArkaitz_bigdata/irondon/R_libraries")
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

# Interval list 
interval_list <- paste0(Data_shared, Sys.getenv("INTERVAL_LIST"))

# Mapped sequences
mapped_directory <- paste0(Base_dir, Sys.getenv("MAPPED_DIRECTORY"))

# Output PON
pon_directory <- paste0(Base_dir,Sys.getenv("PON_DIRECTORY"))

# Jobs dir
jobs_dir <- paste0(Base_dir, Sys.getenv("JOBS_DIR"))

# Jobs requirements
cpu <- Sys.getenv("CPU_MUTECT2")
memory <- Sys.getenv("MEMORY_MUTECT2")
time <- Sys.getenv("TIME_MUTECT2")
partition <- Sys.getenv("PARTITION_MUTECT2")
mail <- Sys.getenv("MAIL")
################################################################################



################################################################################
###############################  JOBS LOOP  ####################################
################################################################################

for(s in 1:length(samples)){

  # Name of the file (regarding the sample variable)
  filename <- paste(jobs_dir , "03_SNVs_step1_1_PON_Mutect2_", samples[s], ".sh", sep = "") 

  # Input file for tumor sample
  input_tumor <- paste(mapped_directory,samples[s],"_mapped.rg.sorted.dup.bqsr.bam", sep ="")
  
  # Tumor ID
  tumor_id <- paste("Sample_", samples[s], sep ="")

  # Output to save VCF files  
  output_vcf <- paste(pon_directory,samples[s],".vcf.gz", sep ="")

  # File to save the variants detected
  cat(
    "#!/bin/bash",
    paste("#SBATCH --job-name=", samples[s], sep = ""),
    paste("#SBATCH -o logs/03_SNVs_step1_1_PON_Mutect2_", samples[s], ".out", sep = ""),
    paste("#SBATCH -e logs/03_SNVs_step1_1_PON_Mutect2_", samples[s], ".err", sep = ""),
    paste("#SBATCH --cpus-per-task=", cpu, sep = ""),
    paste("#SBATCH --mem=", memory, sep = ""),
    paste("#SBATCH --time=", time, sep = ""),
    paste("#SBATCH --partition=", partition, sep = ""),
    paste("#SBATCH --mail-user=",mail,sep = ""),
    "#SBATCH --mail-type=FAIL",

    "\n\n",
    paste("source ",conda_activate, sep = ""), 
    paste("conda activate ", conda_env_wes,sep = ""), 
    "\n\n",

    paste("gatk Mutect2 -R", genome, 
    "-I", input_tumor,
    "-tumor", tumor_id, 
    "--germline-resource", snps_indels_data,
    "-L", interval_list,
    "-O", output_vcf, sep = " "),

    file = filename, sep = "\n", append = FALSE)

  system(paste("sbatch", filename, sep = " "))

}