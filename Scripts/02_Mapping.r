#/opt/R/R-4.1.2/bin/R

################################################################################
####### THIS CODES CREATES A SCRIPT FOR EACH FILE FASTA AND SUMIT JOBS #########
################################################################################

#------------------------------ ALIGNMENT WITH BWA -----------------------------
# BWA (Burrows-Wheeler Aligner) is a popular tool for performing sequence alignment.
# Before running the alignment tool, it is escential to index the reference genome:

# bwa index reference.fasta
# Indexing the reference genome is a critical preprocessing step that enhances the 
# efficiency, speed, and accuracy of the sequence alignment process, making it essential 
# for effective analysis of next-generation sequencing data.

# ADD READ GROUPS BEFORE MARKING DUPLICATES
# Adding read group information to BAM files enhances the accuracy, reliability, and interpretability 
# of sequencing data analysis. 
#-------------------------------------------------------------------------------

#-------------------------------- MARK DUPLICATES ------------------------------
# Why Mark Duplicates?
# 1. Removal of PCR Duplicates
# 2. Avoidance of Bias
# 3. Improvement of Variant Calling
# When to Remove Duplicates?
# While marking duplicates is important for quality control and accurate analysis, whether to remove 
# them entirely depends on the specific analysis and the goals of the study:
# 1. Variant Calling
# 2. Coverage Analysis: If the goal of your analysis is to assess coverage depth or other metrics that 
# may be influenced by duplicate reads, you may choose to retain duplicates for these analyses.
# 3. Downstream Analyses: Certain downstream analyses, such as allele-specific expression or haplotype 
# phasing, may require the retention of duplicates to accurately capture the underlying biology.
# 4. Specific Requirements: Some applications or guidelines may require the removal of duplicates for 
# specific analyses or to meet certain quality standards. In such cases, it's important to carefully 
# consider the implications of removing duplicates and ensure that it aligns with the goals of the study.
#-------------------------------------------------------------------------------

#-------------------------- BASE RECALIBRATION ALGORITHM ------------------------

# The last step of pre-processing mapped reads is the base quality score recalibration 
# (BQSR) stage. 

# Why base recalibration is necessary?

# 1. Correction of Systematic Errors
# 2. Improved Variant Calling
# 3. Enhanced Mapping Quality
# 4. Consistency Across Samples
# 5. Quality Control
# BQSR process applies machine learning to model the possible errors and adjust the base 
# quality scores accordingly.

# A step-by-step guide on how to perform base calibration with GATK4:

# 1. Prepare Input Files: Before running base recalibration, make sure you have the following input files:

#   * Aligned BAM file (e.g., aligned_reads.bam)
#   * Known sites of variation (e.g., dbSNP VCF file, dbsnp.vcf)
#   * Reference genome FASTA file (e.g., reference.fasta)

# 2. Index the Compressed Known Sites
# If your VCF files are compressed as .vcf.gz, make sure they are indexed using tabix:

# tabix -p vcf known_sites.vcf.gz

# This step is crucial because GATK requires indexed VCF files to quickly access the data.

# 3. Run BaseRecalibrator: Use the BaseRecalibrator tool to build a recalibration 
# model based on the aligned BAM file and known sites of variation.

# 4. Apply Base Quality Score Recalibration (BQSR): Use the ApplyBQSR tool to apply the recalibration 
# model to the original BAM file and generate a recalibrated BAM file.
#-------------------------------------------------------------------------------

################################################################################


################################################################################
###########################  DATA CONFIGURATION  ###############################
################################################################################

# Configuration env file
conf_env <- "../00_conf_env.env"

# Loading necessary libraries and data
.libPaths("/vols/GPArkaitz_bigdata/irondon/R_libraries")
library("dotenv")
load_dot_env(file = conf_env)  # Loading variables from .env file

# Setting based list directories
Base_dir <- Sys.getenv("BASE_DIR")
Data_shared <- Sys.getenv("DATA_SHARED")

# FASTQs directories
FASTQs_directory <- paste0(Data_shared, Sys.getenv("FASTQs_DIR"))

# List with the names of the files we will process, only keeping the sample name.
samples <- list.files(path = FASTQs_directory, pattern = "*_trimmed_1.fastq")
samples <- gsub("_trimmed_1.fastq.gz", "", samples)

# Conda Activate
conda_activate <- Sys.getenv("CONDA_ACTIVATE")

# Conda environment 
conda_env_wes <- paste0(Data_shared, Sys.getenv("CONDA_ENV_WES"))

# Mouse genome
genome <- paste0(Data_shared, Sys.getenv("GENOME"))

# Germline source
snps_indels_data <- paste0(Data_shared,Sys.getenv("SNPS_INDELS_DATA"))

# Interval list 
interval_list <- paste0(Data_shared, Sys.getenv("INTERVAL_LIST"))

# Bed file
Bed_file <-  paste0(Data_shared,Sys.getenv("BED_FILE"))

# Mapped sequences directory
Mapped_directory <-  paste0(Base_dir, Sys.getenv("MAPPED_DIRECTORY"))

# Duplicate metrics directory
Marked_dup_Directory <- paste0(Base_dir,Sys.getenv("MARKED_DUPLICATES_DIR"))

# Recall directory
Recall_directory <-  paste0(Base_dir,Sys.getenv("RECALL_DIR"))

# Jobs dir
jobs_dir <- paste0(Base_dir, Sys.getenv("JOBS_DIR"))

# Mapping QC directory
mapping_qc_dir <- paste0(Base_dir, Sys.getenv("MAPPED_QC_DIRECTORY"))

# Coverage file
Coverage_file <- paste0(Base_dir, Sys.getenv("COVERAGE_FILE"))

# Jobs requirements
cpu <- Sys.getenv("CPU_MAPPING")
memory <- Sys.getenv("MEMORY_MAPPING")
time <- Sys.getenv("TIME_MAPPING")
partition <- Sys.getenv("PARTITION_MAPPING")
mail <- Sys.getenv("MAIL")

################################################################################


################################################################################
###############################  JOBS LOOP  ####################################
################################################################################

for(s in 1:length(samples)){

  # Name of the file (regarding the sample variable)
  filename <- paste(jobs_dir , "02_Mapping_", samples[s], ".sh", sep = "") 

  # Forward strand
  input_fastq_read1 <- paste(FASTQs_directory, samples[s], "_trimmed_1.fastq.gz", sep = "")

  # Reverse strand
  input_fastq_read2 <- paste(FASTQs_directory, samples[s], "_trimmed_2.fastq.gz", sep = "")
  
  # Output bam
  output_bam <- paste(Mapped_directory, samples[s], "_mapped.bam", sep = "") 

  # Output directory to store with added read groups
  output_rg_bam <- paste(Mapped_directory, samples[s], "_mapped.rg.bam", sep = "")

  # Output directory to store sorted bam files
  output_rg_sorted_bam <- paste(Mapped_directory, samples[s], "_mapped.rg.sorted.bam", sep = "") 

  # Output directory to store sorted bam files
  output_rg_sorted_dup_bam <- paste(Mapped_directory, samples[s], "_mapped.rg.sorted.dup.bam", sep = "") 

  # Output directory to store sorted bam files
  output_rg_sorted_dup_bqsr_bam <- paste(Mapped_directory, samples[s], "_mapped.rg.sorted.dup.bqsr.bam", sep = "") 

  # Recal datafile
  recal_data_file <- paste(Recall_directory,"recal_data_",samples[s],".table", sep ="")

  cat(
    "#!/bin/bash",
    paste("#SBATCH --job-name=", samples[s], sep = ""),
    paste("#SBATCH -o logs/02_Mapping_", samples[s], ".out", sep = ""),
    paste("#SBATCH -e logs/02_Mapping_", samples[s], ".err", sep = ""),
    paste("#SBATCH --cpus-per-task=", cpu, sep = ""),
    paste("#SBATCH --mem=", memory, sep = ""),
    paste("#SBATCH --time=", time, sep = ""),
    paste("#SBATCH --partition=", partition, sep = ""),
    paste("#SBATCH --mail-user=",mail,sep = ""),
    "#SBATCH --mail-type=FAIL",

    "\n\n",
    paste("source ",conda_activate, sep = ""), 
    paste("conda activate ", conda_env_wes, sep = ""), 
    "\n\n",

    #"#---------------------------------- ALIGNMENT WITH BWA ----------------------------------",
    "\n",
    "# STEP 1 and 2: BWA TO ALIGN FASTQ FILES WITHOUT SAVING SAM FILES",
    paste("bwa mem -t ",cpu, " ", genome," ", input_fastq_read1, " ", input_fastq_read2, " | samtools view -@",cpu," -S -b -o ", output_bam, sep = ""),

    "\n",
    "# STEP 3: Add read groups",
    paste("picard AddOrReplaceReadGroups I=", output_bam," O=", output_rg_bam, " RGID=Sample_",samples[s], " RGLB=Library1 RGPL=illumina RGPU=Unit1 RGSM=Sample_",samples[s], sep = ""),

    "\n",
    "# STEP 4: Sort BAM USING SAMTOOLS",
    paste("samtools sort -@",cpu," ", output_rg_bam, " -o ", output_rg_sorted_bam, sep = ""),
    paste("samtools index", output_rg_sorted_bam, sep = " "),

    "\n\n",
   
    "#------------------------------------ MARK DUPLICATES -----------------------------------",
    "\n",

    paste("picard MarkDuplicates INPUT=", output_rg_sorted_bam," OUTPUT=", output_rg_sorted_dup_bam, " METRICS_FILE=", Marked_dup_Directory , "marked_dup_metrics_",samples[s], ".txt CREATE_INDEX=TRUE REMOVE_DUPLICATES=TRUE", sep = ""),

    "\n\n",

    "#--------------------------- BASE QUALITY SCORE RECALIBRATION ---------------------------",
    "\n",

    "# Step 1: Build the model",
    paste("gatk BaseRecalibrator -I", output_rg_sorted_dup_bam, 
    "-R", genome, 
    "-L", interval_list, 
    "--known-sites", snps_indels_data, 
    "-O", recal_data_file, sep = " "),
    "\n",

    "# Step 2: Apply the model to adjust the base quality scores",
    paste("gatk ApplyBQSR -I", output_rg_sorted_dup_bam, 
    "-R", genome, 
    "-L", interval_list, 
    "--bqsr-recal-file", recal_data_file,
    "-O", output_rg_sorted_dup_bqsr_bam, sep = " "),
    "\n",
    "# Step 3: Indexing the final results",
    paste("samtools index", output_rg_sorted_dup_bqsr_bam, sep = " "),

    "\n\n",

    "#--------------------------------- MAPPING QC METRICS ----------------------------------",
    "\n",
    
    "# Step1: Mapping QC using qualimap",
    paste("mkdir ", mapping_qc_dir,"qualimap_",samples[s], sep = ""),
    
    paste("samtools flagstat ", output_rg_sorted_dup_bqsr_bam, " > ", mapping_qc_dir,samples[s],"_samtools_qc.txt", sep =""),

    #paste("qualimap bamqc --paired --skip-duplicated --sorted -bam ", output_rg_sorted_dup_bqsr_bam, " -outdir ", mapping_qc_dir,"qualimap_",samples[s],"/", sep = ""),
    
    #"\n",
    #paste("mean_depth=$(bedtools coverage -a ",Bed_file," -b ", output_rg_sorted_dup_bqsr_bam, " -d | awk '{sum+=$5} END {if (NR>0) print sum/NR; else print '0'}')", sep = ""),
    #paste("echo -e", samples[s],"\t${mean_depth} >> ",samples[s], "_",Coverage_file, sep = ""),

    file = filename, sep = "\n", append = FALSE)

  system(paste("sbatch", filename, sep = " "))

}