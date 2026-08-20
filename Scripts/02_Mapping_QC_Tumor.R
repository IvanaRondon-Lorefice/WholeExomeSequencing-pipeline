################################################################################
# MAPPED QUALITY CONTROL AC-82_WESmouse
################################################################################

#| Statistical analysis of flagstat obtained from samtools. An example of an output file
#| is the following:

#| 93472189 + 0 in total (QC-passed reads + QC-failed reads)
#| 92934426 + 0 primary
#| 0 + 0 secondary
#| 537763 + 0 supplementary
#| 0 + 0 duplicates
#| 0 + 0 primary duplicates
#| 93364671 + 0 mapped (99.88% : N/A)
#| 92826908 + 0 primary mapped (99.88% : N/A)
#| 92934426 + 0 paired in sequencing
#| 46449125 + 0 read1
#| 46485301 + 0 read2
#| 92158588 + 0 properly paired (99.17% : N/A)
#| 92793352 + 0 with itself and mate mapped
#| 33556 + 0 singletons (0.04% : N/A)
#| 484746 + 0 with mate mapped to a different chr
#| 299818 + 0 with mate mapped to a different chr (mapQ>=5)

#| Where:
#| Total reads: 93,472,189 reads in total, which includes both QC-passed and QC-failed reads.

#| Primary: 92,934,426 primary alignments. These are the main alignments and do not include 
#| secondary or supplementary alignments. In practical terms, primary alignments are typically 
#| the most reliable and accurate alignments for each read. They are usually chosen from among 
#| multiple possible alignments for a given read, such as when a read can map equally well to 
#| multiple locations in the reference genome.

#| Supplementary: 537,763 supplementary alignments. Supplementary alignments represent chimeric 
#| alignments, such as those spanning fusion junctions.

#| Duplicates: 0 duplicate reads. Duplicate reads are identical reads that likely resulted 
#| from PCR amplification during library preparation.

#| Primary duplicates: 0 primary duplicate reads. These are duplicate reads among the primary alignments.

#| Mapped reads: 93,364,671 reads mapped to the reference genome, representing 99.88% of the total reads.

#| Primary mapped reads: 92,826,908 primary reads mapped to the reference genome, also 
#| representing 99.88% of the total reads.

#| Paired in sequencing: 92,934,426 reads were part of paired-end sequencing.

#| Read1: 46,449,125 reads are from the first mate (read1).
#| Read2: 46,485,301 reads are from the second mate (read2).

#| Properly paired: 92,158,588 reads are properly paired, representing 99.17% of the paired reads.

#| With itself and mate mapped: 92,793,352 reads have both themselves and their mates mapped to the reference genome.

#| Singletons: 33,556 reads are singletons, meaning they are not part of a properly paired alignment.

#| With mate mapped to a different chr: 484,746 reads have their mates mapped to a different chromosome.

#| With mate mapped to a different chr (mapQ>=5): 299,818 reads have their mates mapped to a 
#| different chromosome with a mapping quality score (mapQ) of at least 5.

#| Secondary alignments refer to additional alignments for a read that aligns to multiple locations 
#| in the genome. These occur when a read maps well to more than one location. The "primary alignment" 
#| is the best-scoring alignment, while "secondary alignments" are alternative locations where the read 
#| could map.
#| Secondary alignment = 0 indicates that none of the reads in these samples had alternative alignments. 
#| In other words, for each read, only one location in the genome was considered the best and reported 
#| as the primary alignment.

################################################################################


################################################################################
# LIBRARIES AND DATA
################################################################################
suppressMessages(library(ggplot2))
suppressMessages(library(dotenv))

#| For plots
theme_set(theme_classic())

#| Setting working dir
setwd("X:/irondon/AC-82_WESmouse_pipeline")

#| Configuration env file
load_dot_env(file = "00_conf_env.env") 

#| List with the names of the files we will process, only keeping the sample name.
samples <- c(strsplit(Sys.getenv("NORMAL_PON_SAMPLES"), ",")[[1]], strsplit(Sys.getenv("TUMOR_PON_SAMPLES"), ",")[[1]])

#| Mapping QC directory
mapping_qc_dir <- paste0(Sys.getenv("R_BASE_DIR"), Sys.getenv("MAPPED_QC_DIRECTORY"))
dup_metrics_dir <- paste0(Sys.getenv("R_BASE_DIR"), Sys.getenv("MARKED_DUPLICATES_DIR"))

#| Condition
cond <- "PON"
################################################################################


################################################################################
# STATISTICS FROM FLAGSTATS
################################################################################
metrics_name <- c("Total", 
                  "Primary", 
                  "Secondary", 
                  "Supplementary",
                  "Duplicates",
                  "Primary duplicates",
                  "Mapped", 
                  "Primary mapped", 
                  "Paired in sequencing", 
                  "Read1", 
                  "Read2", 
                  "Properly Paired", 
                  "Mate mapped",
                  "Singletons", 
                  "Different Chr", 
                  "Different Chr mapQ>=5")

#| Creating dataframe to save metrics per samples
metrics_df <- data.frame(row.names = metrics_name)

#| Loop over all samples
for (sample in samples){
  
  flagstat_file <- readLines(paste(mapping_qc_dir, sample,"_samtools_qc.txt", sep = ""))
  
  metrics_df[[sample]] <- gsub(" +.*", "", flagstat_file)
}
################################################################################


################################################################################
# STATISTICS FROM DUPLICATES
################################################################################

duplicates <- readLines(paste0(dup_metrics_dir,"marked_dup_metrics.txt"))
col_names <- duplicates[7]
dup_dataframe <- data.frame(row.names = unlist(strsplit(col_names, split = "\t")))

#| Loop over all samples
for (sample in samples){
  
  duplicates <- readLines(paste0(dup_metrics_dir,"marked_dup_metrics_",sample,".txt", sep =""))
  
  dup_dataframe[[sample]] <- unlist(strsplit(duplicates[8], split = "\t"))
}

################################################################################


################################################################################
# PLOTTING RESULTS BY SAMPLE
################################################################################

Samples <- colnames(metrics_df)

#| Computing stats

Primary_mean <- (as.numeric(metrics_df[2,])/as.numeric(metrics_df[1,]))*100
df <- data.frame(value = Primary_mean,
                 Samples = Samples)
ggplot(df, aes(x = reorder(Samples, -value), y =value, fill = Samples))+
  geom_bar(stat="identity", color ="black", fill ="darkorange") +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(angle = 60, hjust=1)) +
  ylab("Primary read (%)")+
  xlab("Samples")
ggsave( paste0(mapping_qc_dir,"Primary_read_",cond,".pdf"), height = 4.5, width = 8.5)

Secondary_mean <- (as.numeric(metrics_df[3,])/as.numeric(metrics_df[1,]))*100
df <- data.frame(value = Secondary_mean,
                 Samples = Samples)
ggplot(df, aes(x = reorder(Samples, -value), y =value, fill = Samples))+
  geom_bar(stat="identity", color ="black", fill ="midnightblue") +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(angle = 60, hjust=1)) +
  ylab("Secondary read (%)")+
  xlab("Samples")
ggsave( paste0(mapping_qc_dir,"Secondary_read_",cond,".pdf"), height = 4.5, width = 8.5)

Supplementary_mean <- (as.numeric(metrics_df[4,])/as.numeric(metrics_df[1,]))*100
df <- data.frame(value = Supplementary_mean,
                 Samples = Samples)
ggplot(df, aes(x = reorder(Samples, -value), y =value, fill = Samples))+
  geom_bar(stat="identity", color ="black", fill ="violet") +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(angle = 60, hjust=1)) +
  ylab("Supplementary (%)")+
  xlab("Samples")
ggsave( paste0(mapping_qc_dir,"Supplementary_percentage_",cond,".pdf"), height = 4.5, width = 8.5)


Supplementary <- as.numeric(metrics_df[4,])
df <- data.frame(value = Supplementary,
                 Samples = Samples)
ggplot(df, aes(x = reorder(Samples, -value), y =value, fill = Samples))+
  geom_bar(stat="identity", color ="black", fill ="violet") +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(angle = 60, hjust=1)) +
  ylab("Supplementary")+
  xlab("Samples")
ggsave( paste0(mapping_qc_dir,"Supplementary_",cond,".pdf"), height = 4.5, width = 8.5)


Mapped_mean <- (as.numeric(metrics_df[7,])/as.numeric(metrics_df[1,]))*100
df <- data.frame(value = Mapped_mean,
                 Samples = Samples)
ggplot(df, aes(x = reorder(Samples, -value), y =value, fill = Samples))+
  geom_bar(stat="identity", color ="black", fill ="green") +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(angle = 60, hjust=1)) +
  ylab("Mapped read (%)")+
  xlab("Samples")
ggsave( paste0(mapping_qc_dir,"Mapped_reads_",cond,".pdf"), height = 4.5, width = 8.5)


Primary_mapped_mean <- (as.numeric(metrics_df[8,])/as.numeric(metrics_df[1,]))*100
df <- data.frame(value = Primary_mapped_mean,
                 Samples = Samples)
ggplot(df, aes(x = reorder(Samples, -value), y =value, fill = Samples))+
  geom_bar(stat="identity", color ="black", fill ="darkgreen") +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(angle = 60, hjust=1)) +
  ylab("Primary mapped reads (%)")+
  xlab("Samples")
ggsave( paste0(mapping_qc_dir,"Primary_mapped_reads_",cond,".pdf"), height = 4.5, width = 8.5)


Paired_in_sequencing_mean <- (as.numeric(metrics_df[9,])/as.numeric(metrics_df[1,]))*100
Properly_paired_mean <- (as.numeric(metrics_df[12,])/as.numeric(metrics_df[1,]))*100
Mate_mapped_mean <- (as.numeric(metrics_df[13,])/as.numeric(metrics_df[1,]))*100
Singlentons_mean <- (as.numeric(metrics_df[14,])/as.numeric(metrics_df[1,]))*100
Diff_Chr_mean <- (as.numeric(metrics_df[15,])/as.numeric(metrics_df[1,]))*100
Diff_Chr_mean_Q <- (as.numeric(metrics_df[16,])/as.numeric(metrics_df[1,]))*100

Duplicates_mean <- (as.numeric(dup_dataframe[9,]))*100
df <- data.frame(value = Duplicates_mean,
                 Samples = Samples)
ggplot(df, aes(x = reorder(Samples, -value), y =value, fill = Samples))+
  geom_bar(stat="identity", color ="black", fill ="red") +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(angle = 60, hjust=1)) +
  ylab("Duplicates (%)")+
  xlab("Samples")
ggsave( paste0(mapping_qc_dir,"Duplicates_",cond,".pdf"), height = 4.5, width = 8.5)



for (i in 1:length(Samples)){
 
  df_plots <- data.frame(metric = c(paste("Primary (", round(Primary_mean[i], 2),"%)", sep =""), 
                                    paste("Secondary (", round(Secondary_mean[i], 2),"%)", sep =""), 
                                    paste("Supplementary (", round(Supplementary_mean[i], 2),"%)", sep =""),
                                    paste("Mapped (", round(Mapped_mean[i], 2),"%)", sep =""), 
                                    paste("Primary mapped (", round(Primary_mapped_mean[i], 2),"%)", sep =""), 
                                    paste("Paired in sequencing (", round(Paired_in_sequencing_mean[i], 2),"%)", sep =""), 
                                    paste("Properly Paired (", round(Properly_paired_mean[i], 2),"%)", sep =""), 
                                    paste("Mate mapped(", round(Mate_mapped_mean[i], 2),"%)", sep ="") ,
                                    paste("Singletons (", round(Singlentons_mean[i], 2),"%)", sep =""), 
                                    paste("Different Chr (", round(Diff_Chr_mean[i], 2),"%)", sep =""), 
                                    paste("Different Chr mapQ>=5 (", round(Diff_Chr_mean_Q[i], 2),"%)", sep =""),
                                    paste("Duplicates (", round(Duplicates_mean[i], 2),"%)", sep ="")),
                         value = c(Primary_mean[i],
                                   Secondary_mean[i],
                                   Supplementary_mean[i],
                                   Mapped_mean[i],
                                   Primary_mapped_mean[i],
                                   Paired_in_sequencing_mean[i],
                                   Properly_paired_mean[i],
                                   Mate_mapped_mean[i],
                                   Singlentons_mean[i],
                                   Diff_Chr_mean[i],
                                   Diff_Chr_mean_Q[i],
                                   Duplicates_mean[i]))
  
  values <- round(df_plots$value[order(df_plots$value, decreasing=T)],2)
  
  #| PLOTTING BY SAMPLE
  ggplot(df_plots, aes(x =value, y = reorder(metric, value))) +
    geom_bar(stat ="identity", fill ="midnightblue")+
    theme(text=element_text(size=16,  family="sans"), 
          legend.key.size = unit(2, 'cm'), 
          plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
    xlab("Percentage (%)") +
    ylab("") 
  ggsave(paste0(mapping_qc_dir,"2_MAPPED_QC_",Samples[i],"_",cond,".pdf"), height = 4.7, width =6.5)
  
  #| SAVING VALUES BY SAMPLE IN EXCEL FILE
  writexl::write_xlsx(df_plots, paste0(mapping_qc_dir,"2_MAPPED_QC_",Samples[i],"_",cond,".xlsx"))
}

################################################################################



################################################################################
# PLOTTING RESULTS ALL
################################################################################

#| Computing stats
Primary_mean <- mean(as.numeric(metrics_df[2,])/as.numeric(metrics_df[1,]))*100
Secondary_mean <- mean(as.numeric(metrics_df[3,])/as.numeric(metrics_df[1,]))*100
Supplementary_mean <- mean(as.numeric(metrics_df[4,])/as.numeric(metrics_df[1,]))*100
Mapped_mean <- mean(as.numeric(metrics_df[7,])/as.numeric(metrics_df[1,]))*100
Primary_mapped_mean <- mean(as.numeric(metrics_df[8,])/as.numeric(metrics_df[1,]))*100
Paired_in_sequencing_mean <- mean(as.numeric(metrics_df[9,])/as.numeric(metrics_df[1,]))*100
Properly_paired_mean <- mean(as.numeric(metrics_df[12,])/as.numeric(metrics_df[1,]))*100
Mate_mapped_mean <- mean(as.numeric(metrics_df[13,])/as.numeric(metrics_df[1,]))*100
Singlentons_mean <- mean(as.numeric(metrics_df[14,])/as.numeric(metrics_df[1,]))*100
Diff_Chr_mean <- mean(as.numeric(metrics_df[15,])/as.numeric(metrics_df[1,]))*100
Diff_Chr_mean_Q <- mean(as.numeric(metrics_df[16,])/as.numeric(metrics_df[1,]))*100
Duplicates_mean <- mean(as.numeric(dup_dataframe[9,]))*100

df_plots <- data.frame(metric = c(paste("Primary (", round(Primary_mean, 2),"%)", sep =""), 
                                  paste("Secondary (", round(Secondary_mean, 2),"%)", sep =""), 
                                  paste("Supplementary (", round(Supplementary_mean, 2),"%)", sep =""),
                                  paste("Mapped (", round(Mapped_mean, 2),"%)", sep =""), 
                                  paste("Primary mapped (", round(Primary_mapped_mean, 2),"%)", sep =""), 
                                  paste("Paired in sequencing (", round(Paired_in_sequencing_mean, 2),"%)", sep =""), 
                                  paste("Properly Paired (", round(Properly_paired_mean, 2),"%)", sep =""), 
                                  paste("Mate mapped(", round(Mate_mapped_mean, 2),"%)", sep ="") ,
                                  paste("Singletons (", round(Singlentons_mean, 2),"%)", sep =""), 
                                  paste("Different Chr (", round(Diff_Chr_mean, 2),"%)", sep =""), 
                                  paste("Different Chr mapQ>=5 (", round(Diff_Chr_mean_Q, 2),"%)", sep =""),
                                  paste("Duplicates (", round(Duplicates_mean, 2),"%)", sep ="")),
                        value = c(Primary_mean,
                                  Secondary_mean,
                                  Supplementary_mean,
                                  Mapped_mean,
                                  Primary_mapped_mean,
                                  Paired_in_sequencing_mean,
                                  Properly_paired_mean,
                                  Mate_mapped_mean,
                                  Singlentons_mean,
                                  Diff_Chr_mean,
                                  Diff_Chr_mean_Q,
                                  Duplicates_mean))

values <- round(df_plots$value[order(df_plots$value, decreasing=T)],2)
ggplot(df_plots, aes(x =value, y = reorder(metric, value))) +
  geom_bar(stat ="identity", fill ="midnightblue")+
  theme(text=element_text(size=16,  family="sans"), 
        legend.key.size = unit(2, 'cm'), 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("Percentage (%)") +
  ylab("") 
ggsave(paste0(mapping_qc_dir,"2_MAPPED_QC_",cond,".pdf"), height = 4.7, width =6.5)
################################################################################

