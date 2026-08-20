#/opt/ohpc/pub/apps/R/R-4.2.1/bin/R


###############################################################################
####### THIS CODES CREATES A SCRIPT FOR EACH FASTQ FILE AND SUMIT JOBS ########
#######   TO PERFORM TRIMMING ON THE FASTQ FILES OF THE SECOND READ    ########
###############################################################################

#| From the cutadapt documentation (paired-end sequence)
#|     cutadapt [options] -o output 1 -p output2 input1 input2
#| [options]: 
#|     -a: Cut adapters of the 3' of Read 1
#|     -A: Cut adapters of the 3' of Read 2
#|     -m: Defines a threshold of the minimum lenght sequence to keep
#|     -q: Defines a threshold of the quality for each sequence
#|     -u: Cut a specific number of nucleotides of Read 1
#|     -U: Cut a specific number of nucleotides of Read 2

###############################################################################


################################################################################
###############################  DATA LIST  ####################################
################################################################################

# Configuration env file
conf_env <- "../00_conf_env.env"

# Loading necessary libraries and data
.libPaths("../R_libraries")
library("dotenv")
load_dot_env(file = conf_env)  # Loading variables from .env file

# Setting based list directories
Base_dir <- Sys.getenv("BASE_DIR")
Data_shared <- Sys.getenv("DATA_SHARED")

# FASTQs directories
FASTQs_directory <- paste0(Data_shared, Sys.getenv("FASTQs_DIR"))

# List with the names of the files we will process, only keeping the sample name.
samples <- list.files(path = paste0(Sys.getenv("DATA_SHARED"),Sys.getenv("PROJECT_NAME"),"/FASTQs/" ), pattern = "*.fastq.gz")
samples <- gsub(".fastq.gz", "", samples)

# FASTQs files directory
FASTQs_directory <- paste0(Sys.getenv("DATA_SHARED"),Sys.getenv("PROJECT_NAME"),"/FASTQs/" )

# FASTQs files trimmed directory
FASTQs_directory_trimmed <- paste0(Data_shared,Sys.getenv("FASTQs_DIR"))

# Creating FASTQs trimmed directory
mkdir(paste0(Data_shared,Sys.getenv("FASTQs_DIR")))

# Conda Activate
conda_activate <- Sys.getenv("CONDA_ACTIVATE")

# Conda environment 
conda_env_wes <- paste0(Data_shared, Sys.getenv("CONDA_ENV_WES"))

# Jobs dir
jobs_dir <- paste0(Base_dir, Sys.getenv("JOBS_DIR"))

# Jobs requirements
cpu <- Sys.getenv("CPU_TRIMMING")
memory <- Sys.getenv("MEMORY_TRIMMING")
time <- Sys.getenv("TIME_TRIMMING")
partition <- Sys.getenv("PARTITION_TRIMMING")
mail <- Sys.getenv("MAIL")
################################################################################


################################################################################
###############################  JOBS LOOP  ####################################
################################################################################
# This creates a loop over all the samples and create a .sh file which is then 
#submitted as a job
for(s in 1:length(samples)){

  # Name of the file (regarding the sample variable)
  filename <- paste(jobs_dir, "TRIMMING_", samples[s], ".sh", sep = "") 

  # Forward strand  (Read 1)
  input1 <- paste(FASTQs_directory, samples[s], "_1.fastq.gz", sep = "")
  input2 <- paste(FASTQs_directory, samples[s], "_2.fastq.gz", sep = "")

  # Output directories
  output1 <- paste(FASTQs_directory_trimmed, samples[s], "_trimmed_1.fastq.gz", sep = "") 
  output2 <- paste(FASTQs_directory_trimmed, samples[s], "_trimmed_2.fastq.gz", sep = "")
  
  cat(
    "#!/bin/sh",
    paste("#SBATCH --job-name=", samples[s], sep = ""),
    paste("#SBATCH -o logs/", samples[s], "_TRIM.out", sep = ""),
    paste("#SBATCH -e logs/", samples[s], "_TRIM.err", sep = ""),
    paste("#SBATCH --cpus-per-task=", cpu, sep = ""),
    paste("#SBATCH --mem=", memory, sep = ""),
    paste("#SBATCH --time=", time, sep = ""),
    paste("#SBATCH --partition=", partition, sep = ""),
    "#SBATCH --mail-user=irondon@cicbiogune.es",
    "#SBATCH --mail-type=FAIL",
    "\n\n",
    paste("source ",conda_activate, sep = ""), 
    paste("conda activate ", conda_env_wes, sep = ""),
    "\n",
    paste("cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -j 0 -q 10 -m 20 --pair-filter=any -o", output1, "-p", output2, input1, input2, sep = " "),
    file = filename, sep = "\n", append = FALSE)
  system(paste("sbatch", filename, sep = " "))

  
}
