#!/bin/bash

#################################################################################
### QUALITY CONTROL ON THE FASTQ FILE WITH FASTQC SOFTWARE FOR RNASEQ ANALYSIS ##
###    READ 1 
#################################################################################

# FastQC reads a set of sequence files and produces from each one a quality control report
#consisting of a number of different modules, each one of which will help to identify a di-
#fferent potential type of problem in your data. Basically, the command line has to follow:
#
#    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
#          [-c contaminant file] seqfile1 .. seqfileN
#where:
#   -h --help:     Print the help file
#   -v --version:  Print the version of the program
#   -o --outdir:   Creates an output folder with the same name where all the outfiles will
#be stored
#   --casava:      Files come from raw casava output. Files in the same sample group (diffe-
#ring only by the group number) will be analysed as a set rather than individually. Files 
#must have the same names given to them by casava (including being gzipped and ending with 
#.gz) otherwise they won't be grouped together correctly.
#   --nano:        Files come from nanopore sequences and are in the fast5 format
#   -t --threads   Specifies the number of files which can be processed simultaneously.
#Take into account that each thread will allocate 250MB of memory, so you shouldn't run 
#more threads than your available memory will cope with, and not more than 6 threads on a 
#32 bit machine

##################################################################################


#SBATCH --job-name="FASTQC_READ_1"
#SBATCH -o logs/FASTQC_READ1.out
#SBATCH -e logs/FASTQC_READ1.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=3GB
#SBATCH --partition=FAST
#SBATCH --time=12:00:00

# Load environment variables from .env file
set -a
. ../conf_env.env
set +a

# Activate conda environment
source "$CONDA_ACTIVATE"
conda activate "$DATA_SHARED/$CONDA_ENV_WES"

# Define paths using exported variables
FASTQs_dir="${DATA_SHARED}${PROJECT_NAME}/FASTQs_trimmed"
FASTQCs_out="${DATA_SHARED}${PROJECT_NAME}/FASTQCs_trimmed"

# Quality control of the sequences from FASTqs files (all the files in the folder)
fastqc $FASTQs_dir/*_1.fastq.gz -o $FASTQCs_out/ 
