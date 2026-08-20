#/opt/R/R-4.1.2/bin/R

################################################################################
####### THIS CODES CREATES A SCRIPT FOR EACH FILE FASTA AND SUMIT JOBS #########
################################################################################

#----------------------------------- MUTECT2 ----------------------------------

#| Variant calling is the process of identifying genetic variations, such as single nucleotide polymorphisms (SNPs), 
#| insertions, deletions, and structural variants, from DNA sequencing data. It is a crucial step in genomic analysis, 
#| providing insights into genetic differences between individuals and populations. Mutect2, developed by the Broad Institute, 
#| is a widely used tool for somatic variant calling in cancer genomics. It is specifically designed to identify somatic 
#| mutations in tumor-normal paired samples with high sensitivity and specificity, making it ideal for detecting low-frequency 
#| variants in complex tumor genomes. Mutect2's sophisticated algorithms and advanced filtering strategies enable accurate 
#| identification of somatic mutations, helping researchers and clinicians unravel the genetic basis of cancer and guide 
#| personalized treatment decisions.

#------------------------------------------------------------------------------

################################################################################


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

# Only tumor samples
samples <- strsplit(Sys.getenv("TUMOR_PON_SAMPLES"), ",")[[1]]

# Conda Activate
conda_activate <- Sys.getenv("CONDA_ACTIVATE")

# Conda environment WES
conda_env_wes <- paste0(Data_shared, Sys.getenv("CONDA_ENV_WES"))

# Conda environment Variance Effect Predictor
conda_env_vep <-  Sys.getenv("CONDA_ENV_VEP")

# Directory where VEP data is stores
vep_data <- Sys.getenv("VEP_DIR")

# Mouse genome
genome <- paste0(Data_shared, Sys.getenv("GENOME"))

# Germline source
snps_indels_data <- paste0(Data_shared,Sys.getenv("SNPS_INDELS_DATA"))

# Interval list 
interval_list <- paste0(Data_shared, Sys.getenv("INTERVAL_LIST"))

# Mapped sequences
mapped_directory <- paste0(Base_dir, Sys.getenv("MAPPED_DIRECTORY"))

# Input of normals obtained from PON
pon_vcf_gz <- paste0(paste0(Base_dir,Sys.getenv("PON_DIRECTORY")), "PON/PON_min-sample-count-1_2.vcf.gz")

# Tumor directory 
tumor_directory <- paste0(Base_dir, Sys.getenv("TUMOR_PON_DIR"))

# Jobs dir
jobs_dir <- paste0(Base_dir, Sys.getenv("JOBS_DIR"))

# Metrics directory 
metrics_dir <- paste0(Base_dir, Sys.getenv("TUMOR_PON_METRICS_DIR"))

# Tumour annotation directory
tumor_annotation_dir <- paste0(Base_dir, Sys.getenv("TUMOR_PON_ANNOTATIONS_DIR"))

# Specie and Version for VEP
specie <- Sys.getenv("SPECIE_VEP")
cache <- Sys.getenv("VERSION_ENSEMBL")

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
  filename <- paste( jobs_dir, "03_SNVs_step2_Tumor-PON_Mutect2_", samples[s], ".sh", sep = "") 

  # Input file for tumor sample
  input_tumor <- paste(mapped_directory,samples[s],"_mapped.rg.sorted.dup.bqsr.bam", sep ="")
  
  # Tumor ID
  tumor_id <- paste("Sample_", samples[s], sep ="")

  # Output to save VCF files  
  output_vcf_vc <- paste(tumor_directory,samples[s],".vcf.gz", sep ="")

  # Output to save VCF files  
  output_vcf_vc_fmc <- paste(tumor_directory,samples[s],"_FMC.vcf.gz", sep ="")
  
  # Output to save VCF files  
  output_vcf_vc_filtered <- paste(tumor_directory,samples[s],"_filtered.vcf.gz", sep ="")

  # File to save the variants detected
  cat(
    "#!/bin/bash",
    paste("#SBATCH --job-name=", samples[s], sep = ""),
    paste("#SBATCH -o logs/03_SNVs_step2_Tumor-PON_Mutect2_", samples[s], ".out", sep = ""),
    paste("#SBATCH -e logs/03_SNVs_step2_Tumor-PON_Mutect2_", samples[s], ".err", sep = ""),
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
    "--panel-of-normals", pon_vcf_gz,
    "-L ", interval_list, 
    "-O", output_vcf_vc, sep = " "),

    "\n",

    paste("gatk FilterMutectCalls -R", genome, 
          "-V", output_vcf_vc,
          "-O", output_vcf_vc_fmc, sep = " "),
    
    "\n",
    paste("bcftools view -f PASS -O z -o", output_vcf_vc_filtered, output_vcf_vc_fmc, sep = " "),
    paste("bcftools index -t", output_vcf_vc_filtered, sep = " "),
    paste("bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/DP\\t%INFO/ECNT\\t%INFO/POP_AF\\t%INFO/P_GERMLINE\\t%INFO/TLOD\\t%FORMAT\\n' ", output_vcf_vc_filtered, 
    " > ", 
    metrics_dir, samples[s], "_filtered.txt", sep = ""),

    "\n",

    "conda deactivate",
    paste("source ",conda_activate, sep = ""), 
    paste("conda activate ", conda_env_vep,sep = ""), 
    paste("vep --cache --offline --appris --biotype --buffer_size 500 --check_existing --distance 5000 --CACHE_VERSION ",cache," --mane --regulatory --show_ref_allele --sift b --species ",specie," --symbol --transcript_version --tsl --uploaded_allele --force_overwrite --dir ", vep_data,  
          " --fasta ", genome, 
          " --input_file ", output_vcf_vc_filtered,
          " --output_file ",tumor_annotation_dir,samples[s],"_vep.txt", sep = ""),

    paste("sed -n '/^#Uploaded_variation/,$p' ",tumor_annotation_dir,samples[s],"_vep.txt", " | sed '1s/^#//' > ", tumor_annotation_dir, samples[s],".txt",sep = ""), 

    file = filename, sep = "\n", append = FALSE)

  system(paste("sbatch", filename, sep = " "))

}