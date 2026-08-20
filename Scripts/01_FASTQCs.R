################################################################################
#                   SUMMARY OF THE FASTQC ANALYSIS
################################################################################

# Once the fastqc software is applied to the fastq files obtained from the RNASeq
#we need to summarize the results obtained on each files (the zip files) regarding
#sequence quality, adapter content, etc..

# This script is written to evaluate the results from the FastQC applied to the
#MG-04_Illumina_totalRNASeq

################################################################################


################################################################################
#|  LIBRARIES AND DATA 
################################################################################
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))

#| Setting working dir
setwd("../irondon/")

# Configuration env file
load_dot_env(file = "00_conf_env.env") 

#| Setting based list directories
Base_dir <- Sys.getenv("R_BASE_DIR")
Data_shared <- Sys.getenv("R_DATA_SHARED")

#| FASTQs files trimmed directory
FASTQs_directory_trimmed <- paste0(Data_shared,Sys.getenv("FASTQs_DIR"))

#| FASTQCs files trimmed directory
FASTQCs_directory_trimmed <- paste0(Data_shared,Sys.getenv("PROJECT_NAME"),"/FASTQCs_trimmed/")

# Creating FASTQCs trimmed directory
mkdir(paste0(Data_shared,Sys.getenv("PROJECT_NAME"),"/FASTQCs_trimmed/"))

# List with the names of the files we will process, only keeping the sample name.
Samples <- list.files(path = FASTQs_directory_trimmed, pattern = "*.fastq.gz")
Samples <- gsub(".fastq.gz", "", Samples)
Samples <- sort(Samples)

#| Empty dataframe to load all data
fastqc_summary <- data.frame()
################################################################################


################################################################################
#|  LOOP  
################################################################################
#| With a loop I go through the Samples character list
for(s in 1:length(Samples)){
  
  # I read the fastqc summary file, but I have to unzip the directory in the process
  base_table <- read.table(unz(paste(FASTQCs_directory_trimmed, Samples[s], "_fastqc.zip", sep = ""), paste(Samples[s], "_fastqc/summary.txt", sep = "")), sep = "\t")
  
  # I only keep the PASS/FAIL/WARN results
  results <- base_table[,1]
  
  # I Add the name of the sample in the first column of the dataframe (row number is "s")
  fastqc_summary[s,1] <- Samples[s]
  
  # With a for loop, I print every PASS/FAIL/WARN data of that sample in the next columns
  for(k in 1:nrow(base_table)){
    fastqc_summary[s, k+1] <- results[k]
  }
}
################################################################################

################################################################################
#|  FORMATING DATAFRAME  
################################################################################
#| Taking the econd column of the base_table that was recollected in the loop. This column has the names of the different studies
col_names <- base_table[,2]

#| Adding the name of the first column to the object
col_names <- c("FastQC Sample", col_names)

#| Applying the colnames to the data frame
colnames(fastqc_summary) <- col_names
################################################################################

################################################################################
#|  PLOTTING  
################################################################################
#| To only plot the studies with the results (without samples), I first have to transpose the data.frame (what will transform it to a matrix, so I retransforme it into a dataframe to edit later) and delete the first column
fastqc_summary_transpose <- as.data.frame(t(fastqc_summary[-1]))

#| Adding a new column with the rownames called "Studies"
fastqc_summary_transpose$Study <- rownames(fastqc_summary_transpose)

#| Now I have to use the melt() function to put the dataframe in long format by the column "Study"
fastqc_summary_molten <- melt(as.data.table(fastqc_summary_transpose), id.vars="Study")

paste0(Base_dir,"1_FASTQCs/FASTQC_Summary_trimmed.pdf")

pdf(file = paste0(Base_dir,"1_FASTQCs/FASTQC_Summary_trimmed.pdf"))
ggplot(fastqc_summary_molten, 
       aes(x = Study, 
           fill = value)) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("PASS" = "#52bd46", "FAIL" = "#e33a14", "WARN" = "#f5b40f")) +
  labs(y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("FASTQC Summary Overview")
dev.off()
################################################################################
