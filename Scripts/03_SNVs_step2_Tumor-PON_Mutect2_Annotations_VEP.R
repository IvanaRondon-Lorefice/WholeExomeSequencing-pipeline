################################################################################
# ANNOTATING SNVs FOUND IN WES DATA OF PROJECT AC-82_WESmouse_pipeline
################################################################################

#| Description

################################################################################




################################################################################
#| LIBRARY AND DATA
################################################################################
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(gprofiler2))
suppressMessages(library(vctrs))
suppressMessages(library(viridis))
suppressMessages(library(readxl))
suppressMessages(library(colorspace))
suppressMessages(library(ggpubr))
suppressMessages(library(cowplot))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(colorspace))
suppressMessages(library(RColorBrewer))
suppressMessages(library(circlize))
suppressMessages(library(devtools))
suppressMessages(library(trackViewer))
suppressMessages(library(dplyr))
suppressMessages(library(forcats))
suppressMessages(library(patchwork))
suppressMessages(library(maftools))
suppressMessages(library(biomaRt))
suppressMessages(library(tidyverse))
suppressMessages(library(dotenv))

#| For plots
theme_set(theme_classic())

#| Setting working dir
setwd("X:/irondon/AC-82_WESmouse_pipeline")

# Configuration env file
load_dot_env(file = "00_conf_env.env") 

#| Dataset from ensembl
genome <- Sys.getenv("DATASET_ENSEMBL")
version <- as.numeric(Sys.getenv("VERSION_ENSEMBL"))

#| Tumor directory
tumour_dir <- paste0(Sys.getenv("R_BASE_DIR"),Sys.getenv("TUMOR_PON_DIR"))

#| Directory Results
dir.results <- paste0(Sys.getenv("R_BASE_DIR"),Sys.getenv("TUMOR_PON_ANNOTATIONS_DIR"), "Results/")

#| Creating necessary folders
dir.create(dir.results)
dir.create(paste0(dir.results,"Violin_plots"))
dir.create(paste0(dir.results,"Bar_plots"))
dir.create(paste0(dir.results,"Oncoplots"))
dir.create(paste0(dir.results,"Histograms"))

#| Reading mouse information (containing gene type info)
ensembl <- useEnsembl(biomart = "genes", dataset = genome, version = version)
annotations <- getBM(attributes = c("gene_biotype", "external_gene_name", "chromosome_name", "start_position", "end_position"), mart = ensembl)
names(annotations) <- c("gene_type","gene_name", "CONTIG", "START", "END" )
annotations$CONTIG <- paste("chr",annotations$CONTIG, sep ="")

#| Custom background for enrichment analysis
custom_bg <- unique(annotations$gene_name[which(annotations$gene_type == "protein_coding")])

#| Sample information
sample_info <- read_xlsx("X:/irondon/AC-82_WESmouse/Info_mouse.xlsx")
sample_info$Samples <- NA
sample_info$Samples <- paste0("Sample_", sample_info$`Tube ID`)
################################################################################


################################################################################
#| CONFIRMING MUTECT2 METRICS 
################################################################################

confirm_unique <- c("GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:PGT:PID:SA_MAP_AF:SA_POST_PROB", "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB")

#| Extracting sample name
Samples <- list.files(path = tumour_dir, pattern = "*_filtered.vcf.gz.tbi")
Samples <- gsub("_filtered.vcf.gz.tbi", "", Samples)

collect_V10 <- c()

for (i in 1:length(Samples)){

  #| Collecting metrics
  sample_metric <- read.table(paste(tumour_dir, "03_SNVs_step2_1_TUMOR_METRICS_BCFTOOLS/",Samples[i], "_filtered.txt",sep =""), sep ="\t")
  
  #| Saving 
  collect_V10 <- c(collect_V10, unique(sample_metric$V10))
}

collect_V10 <- unique(collect_V10)
collect_V10 <- collect_V10[order(collect_V10)]

if (length(unique(collect_V10)) > 1) {
  
  if (all(confirm_unique == unique(collect_V10))) print("Confirmed values of Mutect2 metrics")
  
} else if (any(confirm_unique == unique(collect_V10)[1])) {
  
  print("Confirmed values of Mutect2 metrics")
  
} else { print("Code will not work") }



################################################################################


################################################################################
#| PROCESSING VEP RESULTS
################################################################################
#| Extracting sample name
Samples <- list.files(path = tumour_dir, pattern = "*_filtered.vcf.gz.tbi")
Samples <- gsub("_filtered.vcf.gz.tbi", "", Samples)

#| Collecting metrics
i <- 1
sample_metric <- read.table(paste(tumour_dir, "03_SNVs_step2_1_TUMOR_METRICS_BCFTOOLS/",Samples[i], "_filtered.txt",sep =""), sep ="\t")

#| Setting the conditions
cond1 <- "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB"
cond2 <- "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:PGT:PID:SA_MAP_AF:SA_POST_PROB"

if( !is.na(unique(sample_metric$V10)[2]) ){
  
  if ( unique(sample_metric$V10)[1] == cond1){
    
    #| Case 1 unique(sample_metric$V10)[1] "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB"
    values_1 <- as.data.frame(do.call(rbind, strsplit(sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[1])], ":")), stringsAsFactors = FALSE)
    names(values_1) <- c("GT", "AD", "AF", "F1R2", "F2R1", "MBQ", "MFRL", "MMQ", "MPOS", "SA_MAP_AF", "SA_POST_PROB")
    values_1$PGT <- NA
    values_1$PID <- NA
    values_1 <- values_1[,order(colnames(values_1))]
    values_1$CHR <- sample_metric$V1[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$POS <- sample_metric$V2[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$REF <- sample_metric$V3[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$ALT <- sample_metric$V4[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$DP <- sample_metric$V5[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$ECNT <- sample_metric$V6[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$POP_AF <- sample_metric$V7[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$P_GERMLINE <- sample_metric$V8[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$TLOD <- sample_metric$V9[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    
    #| Case 2 unique(sample_metric$V10)[2] "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:PGT:PID:SA_MAP_AF:SA_POST_PROB"
    values_2 <- as.data.frame(do.call(rbind, strsplit(sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[2])], ":")), stringsAsFactors = FALSE)
    names(values_2) <- c("GT", "AD", "AF", "F1R2", "F2R1", "MBQ", "MFRL", "MMQ", "MPOS", "PGT", "PID", "SA_MAP_AF", "SA_POST_PROB")
    values_2 <- values_2[,order(colnames(values_2))]
    values_2$CHR <- sample_metric$V1[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$POS <- sample_metric$V2[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$REF <- sample_metric$V3[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$ALT <- sample_metric$V4[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$DP <- sample_metric$V5[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$ECNT <- sample_metric$V6[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$POP_AF <- sample_metric$V7[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$P_GERMLINE <- sample_metric$V8[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$TLOD <- sample_metric$V9[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    
    
  }else{
    
    #| Case 1 unique(sample_metric$V10)[1] "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB"
    values_1 <- as.data.frame(do.call(rbind, strsplit(sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[1])], ":")), stringsAsFactors = FALSE)
    names(values_1) <- c("GT", "AD", "AF", "F1R2", "F2R1", "MBQ", "MFRL", "MMQ", "MPOS", "PGT", "PID", "SA_MAP_AF", "SA_POST_PROB")
    values_1 <- values_1[,order(colnames(values_1))]
    values_1$CHR <- sample_metric$V1[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$POS <- sample_metric$V2[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$REF <- sample_metric$V3[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$ALT <- sample_metric$V4[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$DP <- sample_metric$V5[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$ECNT <- sample_metric$V6[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$POP_AF <- sample_metric$V7[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$P_GERMLINE <- sample_metric$V8[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$TLOD <- sample_metric$V9[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    
    #| Case 2 unique(sample_metric$V10)[2] "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:PGT:PID:SA_MAP_AF:SA_POST_PROB"
    values_2 <- as.data.frame(do.call(rbind, strsplit(sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[2])], ":")), stringsAsFactors = FALSE)
    names(values_2) <- c("GT", "AD", "AF", "F1R2", "F2R1", "MBQ", "MFRL", "MMQ", "MPOS", "SA_MAP_AF", "SA_POST_PROB") 
    values_2$PGT <- NA
    values_2$PID <- NA
    values_2 <- values_2[,order(colnames(values_2))]
    values_2$CHR <- sample_metric$V1[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$POS <- sample_metric$V2[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$REF <- sample_metric$V3[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$ALT <- sample_metric$V4[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$DP <- sample_metric$V5[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$ECNT <- sample_metric$V6[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$POP_AF <- sample_metric$V7[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$P_GERMLINE <- sample_metric$V8[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    values_2$TLOD <- sample_metric$V9[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
    
  }
  
  #| Joining dataframes
  values <- rbind(values_1, values_2)
  values$AF <- as.numeric(values$AF)
  values$MBQ <- as.numeric(values$MBQ)
  values$MMQ <- as.numeric(values$MMQ)
  values$MPOS <- as.numeric(values$MPOS)
  
  #| Separating data: UNIQUE ALLELES
  values_unique <- values[which(!(duplicated(values$POS) | duplicated(values$POS, fromLast = TRUE) )), ]
  values_unique <- values_unique[,c("CHR","POS","REF", "ALT","AD", "AF", "DP", "MMQ", "MPOS","P_GERMLINE", "TLOD")]
  
}else{
  
  #| Case 1 unique(sample_metric$V10)[1] "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB"
  values_1 <- as.data.frame(do.call(rbind, strsplit(sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[1])], ":")), stringsAsFactors = FALSE)
  names(values_1) <- c("GT", "AD", "AF", "F1R2", "F2R1", "MBQ", "MFRL", "MMQ", "MPOS", "SA_MAP_AF", "SA_POST_PROB")
  values_1$PGT <- NA
  values_1$PID <- NA
  values_1 <- values_1[,order(colnames(values_1))]
  values_1$CHR <- sample_metric$V1[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
  values_1$POS <- sample_metric$V2[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
  values_1$REF <- sample_metric$V3[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
  values_1$ALT <- sample_metric$V4[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
  values_1$DP <- sample_metric$V5[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
  values_1$ECNT <- sample_metric$V6[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
  values_1$POP_AF <- sample_metric$V7[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
  values_1$P_GERMLINE <- sample_metric$V8[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
  values_1$TLOD <- sample_metric$V9[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
  
  values <- values_1
  values$AF <- as.numeric(values$AF)
  values$MBQ <- as.numeric(values$MBQ)
  values$MMQ <- as.numeric(values$MMQ)
  values$MPOS <- as.numeric(values$MPOS)
  
  #| Separating data: UNIQUE ALLELES
  values_unique <- values[which(!(duplicated(values$POS) | duplicated(values$POS, fromLast = TRUE) )), ]
  values_unique <- values_unique[,c("CHR","POS","REF", "ALT","AD", "AF", "DP", "MMQ", "MPOS","P_GERMLINE", "TLOD")]
  
}

values_unique$Samples <- Samples[i]
values_unique_list <- list()

#| Adding functional annotation
ensembl_results <- read.table(paste(tumour_dir,"/03_SNVs_step2_2_TUMOR_ANNOTATIONS_VEP/", Samples[i],".txt", sep =""), head =T)
ensembl_results <- ensembl_results %>%
  mutate(Extra = str_split(Extra, ";")) %>%    # Split by semicolon
  unnest(Extra) %>%                            # Expand into multiple rows
  separate(Extra, into = c("key", "value"), sep = "=") %>%  # Separate key and value
  pivot_wider(names_from = key, values_from = value)  # Reshape into wide format
ensembl_results <- ensembl_results[,c("Location", "Consequence", "IMPACT", "SYMBOL", "Gene", "Existing_variation", "SIFT", "BIOTYPE")]
names(ensembl_results) <- c("Location", "Consequence", "IMPACT", "SYMBOL", "Gene", "Existing_variation", "SIFT", "BIOTYPE")
ensembl_results_unique <- ensembl_results
ensembl_results_unique$POS <- do.call(rbind,strsplit(do.call(rbind, strsplit(ensembl_results_unique$Location, "-"))[,1],":"))[,2]

values_unique <- merge(values_unique, ensembl_results_unique, by = "POS")
values_unique_list[[Samples[i]]] <- values_unique

for (i in 2:length(Samples)){
  
  sample_metric <- read.table(paste(tumour_dir,"03_SNVs_step2_1_TUMOR_METRICS_BCFTOOLS/",Samples[i], "_filtered.txt",sep =""), sep ="\t")
  
  if( !is.na(unique(sample_metric$V10)[2]) ){
    
    if ( unique(sample_metric$V10)[1] == cond1){
      
      #| Case 1 unique(sample_metric$V10)[1] "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB"
      values_1 <- as.data.frame(do.call(rbind, strsplit(sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[1])], ":")), stringsAsFactors = FALSE)
      names(values_1) <- c("GT", "AD", "AF", "F1R2", "F2R1", "MBQ", "MFRL", "MMQ", "MPOS", "SA_MAP_AF", "SA_POST_PROB")
      values_1$PGT <- NA
      values_1$PID <- NA
      values_1$CHR <- sample_metric$V1[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
      values_1$POS <- sample_metric$V2[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
      values_1$REF <- sample_metric$V3[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
      values_1$ALT <- sample_metric$V4[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
      values_1$DP <- sample_metric$V5[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
      values_1$ECNT <- sample_metric$V6[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
      values_1$POP_AF <- sample_metric$V7[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
      values_1$P_GERMLINE <- sample_metric$V8[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
      values_1$TLOD <- sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
      
      #| Case 2 unique(sample_metric$V10)[2] "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:PGT:PID:SA_MAP_AF:SA_POST_PROB"
      values_2 <- as.data.frame(do.call(rbind, strsplit(sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[2])], ":")), stringsAsFactors = FALSE)
      names(values_2) <- c("GT", "AD", "AF", "F1R2", "F2R1", "MBQ", "MFRL", "MMQ", "MPOS", "PGT", "PID", "SA_MAP_AF", "SA_POST_PROB")
      values_2 <- values_2[,order(colnames(values_2))]
      values_2$CHR <- sample_metric$V1[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$POS <- sample_metric$V2[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$REF <- sample_metric$V3[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$ALT <- sample_metric$V4[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$DP <- sample_metric$V5[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$ECNT <- sample_metric$V6[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$POP_AF <- sample_metric$V7[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$P_GERMLINE <- sample_metric$V8[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$TLOD <- sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      
      
    }else{
      
      #| Case 1 unique(sample_metric$V10)[1] "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB"
      values_1 <- as.data.frame(do.call(rbind, strsplit(sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[1])], ":")), stringsAsFactors = FALSE)
      names(values_1) <- c("GT", "AD", "AF", "F1R2", "F2R1", "MBQ", "MFRL", "MMQ", "MPOS", "PGT", "PID", "SA_MAP_AF", "SA_POST_PROB")
      values_1 <- values_1[,order(colnames(values_1))]
      values_1$CHR <- sample_metric$V1[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
      values_1$POS <- sample_metric$V2[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
      values_1$REF <- sample_metric$V3[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
      values_1$ALT <- sample_metric$V4[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
      values_1$DP <- sample_metric$V5[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
      values_1$ECNT <- sample_metric$V6[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
      values_1$POP_AF <- sample_metric$V7[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
      values_1$P_GERMLINE <- sample_metric$V8[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
      values_1$TLOD <- sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
      
      
      #| Case 2 unique(sample_metric$V10)[2] "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:PGT:PID:SA_MAP_AF:SA_POST_PROB"
      values_2 <- as.data.frame(do.call(rbind, strsplit(sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[2])], ":")), stringsAsFactors = FALSE)
      names(values_2) <- c("GT", "AD", "AF", "F1R2", "F2R1", "MBQ", "MFRL", "MMQ", "MPOS", "SA_MAP_AF", "SA_POST_PROB") 
      values_2$PGT <- NA
      values_2$PID <- NA
      values_2 <- values_2[,order(colnames(values_2))]
      values_2$CHR <- sample_metric$V1[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$POS <- sample_metric$V2[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$REF <- sample_metric$V3[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$ALT <- sample_metric$V4[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$DP <- sample_metric$V5[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$ECNT <- sample_metric$V6[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$POP_AF <- sample_metric$V7[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$P_GERMLINE <- sample_metric$V8[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      values_2$TLOD <- sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[2])]
      
    }
    
    #| Joining dataframes
    values <- rbind(values_1, values_2)
    values$AF <- as.numeric(values$AF)
    values$MBQ <- as.numeric(values$MBQ)
    values$MMQ <- as.numeric(values$MMQ)
    values$MPOS <- as.numeric(values$MPOS)
    
    #| Separating data: UNIQUE ALLELES
    values_unique <- values[which(!(duplicated(values$POS) | duplicated(values$POS, fromLast = TRUE) )), ]
    values_unique <- values_unique[,c("CHR","POS","REF", "ALT","AD", "AF", "DP", "MMQ", "MPOS","P_GERMLINE", "TLOD")]
    
  }else{
    
    #| Case 1 unique(sample_metric$V10)[1] "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB"
    values_1 <- as.data.frame(do.call(rbind, strsplit(sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10)[1])], ":")), stringsAsFactors = FALSE)
    names(values_1) <- c("GT", "AD", "AF", "F1R2", "F2R1", "MBQ", "MFRL", "MMQ", "MPOS", "SA_MAP_AF", "SA_POST_PROB")
    values_1$PGT <- NA
    values_1$PID <- NA
    values_1 <- values_1[,order(colnames(values_1))]
    values_1$CHR <- sample_metric$V1[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
    values_1$POS <- sample_metric$V2[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
    values_1$REF <- sample_metric$V3[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
    values_1$ALT <- sample_metric$V4[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
    values_1$DP <- sample_metric$V5[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
    values_1$ECNT <- sample_metric$V6[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
    values_1$POP_AF <- sample_metric$V7[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
    values_2$P_GERMLINE <- sample_metric$V8[which( sample_metric$V10 == unique(sample_metric$V10)[1])]
    values_1$TLOD <- sample_metric$V11[which( sample_metric$V10 == unique(sample_metric$V10 )[1])]
    
    values <- values_1
    values$AF <- as.numeric(values$AF)
    values$MBQ <- as.numeric(values$MBQ)
    values$MMQ <- as.numeric(values$MMQ)
    values$MPOS <- as.numeric(values$MPOS)
    
    #| Separating data: UNIQUE ALLELES
    values_unique <- values[which(!(duplicated(values$POS) | duplicated(values$POS, fromLast = TRUE) )), ]
    values_unique <- values_unique[,c("CHR","POS","REF", "ALT","AD", "AF", "DP", "MMQ", "MPOS","P_GERMLINE", "TLOD")]
    
  }
  
  values_unique$Samples <- Samples[i]
  
  #| Adding functional annotation
  ensembl_results <- read.table(paste(tumour_dir, "03_SNVs_step2_2_TUMOR_ANNOTATIONS_VEP/", Samples[i],".txt", sep =""), head =T)
  ensembl_results <- ensembl_results %>%
    mutate(Extra = str_split(Extra, ";")) %>%    # Split by semicolon
    unnest(Extra) %>%                            # Expand into multiple rows
    separate(Extra, into = c("key", "value"), sep = "=") %>%  # Separate key and value
    pivot_wider(names_from = key, values_from = value)  # Reshape into wide format
  ensembl_results <- ensembl_results[,c("Location", "Consequence", "IMPACT", "SYMBOL", "Gene", "Existing_variation", "SIFT", "BIOTYPE")]
  names(ensembl_results) <- c("Location", "Consequence", "IMPACT", "SYMBOL", "Gene", "Existing_variation", "SIFT", "BIOTYPE")
  ensembl_results_unique <- ensembl_results
  ensembl_results_unique$POS <- do.call(rbind,strsplit(do.call(rbind, strsplit(ensembl_results_unique$Location, "-"))[,1],":"))[,2]
  
  values_unique <- merge(values_unique, ensembl_results_unique, by = "POS")
  values_unique_list[[Samples[i]]] <- values_unique
  
}

#| Concatenate list in a dataframe
i <- 1
values_unique_df <- values_unique_list[[Samples[i]]]

for (i in 2:length(Samples)){
  
  values_unique_df <- rbind(values_unique_df,  values_unique_list[[Samples[i]]])
  
}

#| Only protein coding
values_unique_df <- values_unique_df[which(values_unique_df$BIOTYPE == "protein_coding"),]

#| Sorting by high impact
values_unique_df <- values_unique_df %>%
  mutate(IMPACT = factor(IMPACT, levels = c("HIGH", "MODERATE", "LOW","MODIFIER"))) %>%
  arrange(IMPACT, desc(AF)) 

#| Define impact ranking
impact_rank <- c( "HIGH" = 3, "MODERATE" = 2,"LOW" = 1, "MODIFIER" = 0)

#| Add numeric ranking to the impact column
values_unique_df <- values_unique_df %>%
  mutate(Impact_Rank = as.numeric(impact_rank[IMPACT]))

#| Filter to retain only the highest impact variant for each Sample-POS combination
#| Keeps only one variant per Sample-POS, even if multiple variants have the same impact level.
#| Removes duplicates while prioritizing the highest impact.
#| If multiple variants have the same impact, it keeps the first one in the dataset.
values_unique_df <- values_unique_df %>%
  group_by(Samples, POS) %>%
  filter(Impact_Rank == max(Impact_Rank)) %>%
  slice(1) %>%  # This ensures only one row is kept if ties exist
  ungroup()

values_unique_df[which(values_unique_df$SYMBOL == "Psg23"),c("Samples", "SYMBOL", "POS", "AF")]

#| If a gene have several mutations, retain de most consequential 
values_unique_df <- values_unique_df %>%
  group_by(Samples, SYMBOL) %>%
  filter(Impact_Rank == max(Impact_Rank)) %>%
  slice(1) %>%  # This ensures only one row is kept if ties exist
  ungroup()

values_unique_df[which(values_unique_df$SYMBOL == "Psg23"),c("Samples", "SYMBOL", "POS")]


#| Adding mouse information
values_unique_df<- merge(values_unique_df,sample_info, by ="Samples" )
values_unique_df$Age <- gsub("M", "", values_unique_df$Age)
values_unique_df$Age <- as.numeric(values_unique_df$Age)

#| Saving results
writexl::write_xlsx(values_unique_df, paste0(dir.results, "values_unique_results.xlsx"))

#| Allele frequencing distribution
p <- ggplot(values_unique_df, aes(x = Samples, y = AF, fill =Samples)) +
  geom_violin(trim=FALSE) +
  stat_summary(fun.data=mean_sdl, mult=1, 
               geom="pointrange", color="black") +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size =12, angle = 60, hjust=1, color ="black"),
        axis.text.y = element_text(size =12, color ="black"),
        legend.position = "none")+
  scale_fill_manual(values = rainbow(21))
means <- aggregate(AF  ~ Samples + DOB, values_unique_df, mean)
p + geom_text(data = means, aes(x = reorder(Samples, DOB), label = DOB), # Adjust y for spacing
              hjust = -1,vjust=-1, size = 3.3,angle = 90 )
ggsave(paste0(dir.results,"Violin_plots/AF_Samples_Protein_Coding_First_Priorization.pdf"), height = 4.5, width = 8)


################################################################################


################################################################################
#| VARIANCE PRIORIZATION 1:
#|      * A mutation have to be different across samples. Appears only once in a Sample.
#|      * Selecting only HIGH impact mutations
################################################################################

values_unique_df <- read_xlsx(paste0(dir.results, "values_unique_results.xlsx"))

values_unique_df <- values_unique_df[which(values_unique_df$DP > 20),]
values_unique_df <- values_unique_df[which(values_unique_df$AF > 0.005),]
values_unique_df <- values_unique_df[which(values_unique_df$MMQ > 30),]
values_unique_df <- values_unique_df[which(values_unique_df$MPOS >20),]
min(values_unique_df$DP)

filtered_mutations <- values_unique_df[which( (values_unique_df$BIOTYPE == "protein_coding") ),]

filtered_mutations <- filtered_mutations[which( (filtered_mutations$Consequence == "missense_variant") | 
                                                  (filtered_mutations$Consequence == "missense_variant,splice_region_variant") | 
                                                  (filtered_mutations$Consequence == "frameshift_variant") | 
                                                  (filtered_mutations$Consequence == "inframe_deletion") | 
                                                  (filtered_mutations$Consequence == "synonymous_variant") | 
                                                  (filtered_mutations$Consequence == "inframe_insertion") |
                                                  (filtered_mutations$Consequence == "stop_gained") |
                                                  (filtered_mutations$Consequence == "stop_lost") |
                                                  (filtered_mutations$Consequence == "protein_altering_variant") |
                                                  (filtered_mutations$Consequence == "start_lost") |
                                                  (filtered_mutations$Consequence == "coding_sequence_variant") |
                                                  (filtered_mutations$Consequence == "splice_acceptor_variant")  ), ] 

#| Filtering mutations appearing with the same position across sample and selecting the more consequential
filtered_mutations <- filtered_mutations %>%
  group_by(POS) %>%
  mutate(sample_count = n_distinct(Samples)) %>%
  dplyr::filter(sample_count <= 2) %>%
  ungroup()

filtered_mutations_df <- filtered_mutations[,c("Samples", "POS", "SYMBOL", "AF", "DP", "MMQ", "MPOS", "Consequence", "IMPACT", "sample_count","DOB")]

filtered_mutations_df <- filtered_mutations_df[which(!is.na(filtered_mutations_df$Samples)),]

filtered_mutations_df[which(filtered_mutations_df$SYMBOL == "Psg23"),]

#| Allele frequency
p <- ggplot(filtered_mutations_df, aes(x = reorder(Samples,AF), y = AF, fill =Samples)) +
  geom_violin(trim=FALSE) +
  stat_summary(fun.data=mean_sdl, mult=1, 
               geom="pointrange", color="black") +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size =12, angle = 60, hjust=1, color ="black"),
        axis.text.y = element_text(size =12, color ="black"),
        legend.position = "none")+
  scale_fill_manual(values = rainbow(21))+
  geom_hline(yintercept = 0.5, color ="black")
means <- aggregate(AF  ~ Samples + DOB, filtered_mutations_df, mean)
p + geom_text(data = means, aes(x = reorder(Samples, DOB), label = DOB), # Adjust y for spacing
              hjust = -1,vjust=-1, size = 3.3,angle = 90 )

ggsave(paste0(dir.results,"Violin_plots/AF_Samples_Protein_Coding_Second_Priorization.pdf"), height = 4.5, width = 8)

filtered_mutations_df <- filtered_mutations_df[which(filtered_mutations_df$AF < 0.35),]


#############| BAR PLOT CONSEQUENCE
filtered_mutations_df_plot <- filtered_mutations_df
filtered_mutations_df_plot$Consequence[which( !( (filtered_mutations_df_plot$Consequence == "inframe_insertion") | 
                                                   (filtered_mutations_df_plot$Consequence == "frameshift_variant") | 
                                                   (filtered_mutations_df_plot$Consequence == "inframe_deletion") | 
                                                   (filtered_mutations_df_plot$Consequence == "missense_variant") | 
                                                   (filtered_mutations_df_plot$Consequence == "synonymous_variant") | 
                                                   (filtered_mutations_df_plot$Consequence == "stop_gained")) )] <- "other"
ggplot(filtered_mutations_df_plot, aes(x = fct_infreq(Samples), fill =Consequence)) +
  geom_bar()+
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size =12, angle = 60, hjust=1, color ="black"),
        axis.text.y = element_text(size =12, color ="black"))+
  scale_fill_manual(values = rainbow(7))+
  xlab("Samples") +
  ylab("Mutations")
ggsave(paste0(dir.results,"Bar_plots/Barplot_Protein-Coding.pdf"), height = 4.5, width = 8)


filtered_mutations_df_plot[which(filtered_mutations_df_plot$SYMBOL == "Kdm6b"),]
filtered_mutations_df_plot[which(filtered_mutations_df_plot$SYMBOL == "Cdk8"),]
unique(filtered_mutations_df_plot$SYMBOL[which(filtered_mutations_df_plot$sample_count ==2)])

#############| ONCOPLOT
#| Transforming the names of the variants
filtered_mutations_df$Consequence_maftools <- NA
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "missense_variant")] <- "Missense_Mutation"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "synonymous_variant")] <- "Silent"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "inframe_deletion")] <- "In_Frame_Del"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "inframe_insertion")] <- "In_Frame_Ins"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "frameshift_variant")] <- "Frame_Shift_Del"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "stop_gained")] <- "Nonsense_Mutation"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "stop_lost")] <- "Stop_Codon_Del"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "start_lost")] <- "Start_Codon_Del"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "splice_region_variant")] <- "Splice_Site"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "splice_donor_variant")] <- "Splice_Site"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "frameshift_variant,splice_region_variant")] <- "Splice_Site"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "stop_gained,frameshift_variant")] <- "Nonsense_Mutation"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "inframe_insertion,stop_retained_variant")] <- "In_Frame_Ins"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "splice_acceptor_variant")] <- "Splice_Site"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "stop_retained_variant")] <- "Other"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "start_lost,inframe_deletion")] <- "In_Frame_Del"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "protein_altering_variant,incomplete_terminal_codon_variant")] <- "Other"
filtered_mutations_df$Consequence_maftools[which(filtered_mutations_df$Consequence == "splice_donor_variant,splice_donor_5th_base_variant")] <- "Splice_Site"

maf_df <- filtered_mutations_df %>%
  mutate(
    Tumor_Sample_Barcode = Samples,        # Map sample names
    Hugo_Symbol = SYMBOL,                  # Map gene names
    Chromosome = "Unknown",                # Placeholder if chromosome not available
    Start_Position = POS,                  # Use POS as start
    End_Position = POS,                    # Use POS as end for single-base mutations
    VAF = AF,
    Variant_Classification = Consequence_maftools,
    Variant_Type = "SNV",                  # Placeholder for mutation type
    Reference_Allele = "-",                # Placeholder if data not available
    Tumor_Seq_Allele2 = "-",               # Placeholder if data not available
    Protein_Change = NA_character_         # Placeholder if not available
  ) %>%
  dplyr::select(
    Tumor_Sample_Barcode, 
    Hugo_Symbol, Chromosome, Start_Position, End_Position,VAF,
    Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele2, Protein_Change
  )

#| DEFINING COLORS
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Set2')
names(vc_cols) = c(
  'Missense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Frame_Shift_Del',
  'Splice_Site',
  'In_Frame_Del',
  'Nonsense_Mutation'
)

maf_data <- read.maf(maf = maf_df)
maf_data_df <- as.data.frame(maf_data@gene.summary)
rownames(maf_data_df) <- maf_data_df$Hugo_Symbol
maf_data_df <- maf_data_df[,-c(1)]

# Create a new workbook
wb <- createWorkbook()

# List of data frames and corresponding sheet names
data_list <- list(
  Gene_summary = as.data.frame(maf_data@gene.summary),
  Variants_per_sample = as.data.frame(maf_data@variants.per.sample),
  Variants_classification_summary = as.data.frame(maf_data@variant.classification.summary),
  Data = as.data.frame(maf_data@data),
  Summary= as.data.frame(maf_data@summary),
  Silent_mutations = as.data.frame(maf_data@maf.silent)
  
)

# Add each data frame to the workbook as a new sheet
for (sheet_name in names(data_list)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, data_list[[sheet_name]])
}

# Save the workbook
saveWorkbook(wb,  paste0(dir.results,"Gene_Mutations_Summary.xlsx"), overwrite = TRUE)

#| INCLUDING ALLELE FREQUENCY
genes <- rownames(maf_data_df)[which(maf_data_df$MutatedSamples >=2)]
maf_data@gene.summary$Hugo_Symbol
aml_genes <- maf_data@gene.summary$Hugo_Symbol[c(1:length(genes))]
aml_genes_vaf = subsetMaf(maf = maf_data, 
                          genes = aml_genes, 
                          fields = "VAF", mafObj = FALSE)[,mean(VAF, na.rm = TRUE), Hugo_Symbol]

colnames(aml_genes_vaf)[2] = "VAF"

pdf(paste0(dir.results,"Oncoplots/Oncoplots_Protein-Coding.pdf"), height = 4, width = 5)
oncoplot(
  maf = maf_data,
  genes = aml_genes,
  leftBarData = aml_genes_vaf,
  leftBarLims = c(0, 1),
  colors = vc_cols,
  bgCol = "grey90",
  borderCol = "grey100",
  showTumorSampleBarcodes = TRUE,
  fontSize  = 0.5
)
dev.off()


#############| ENRICHMENT
#|      * ALL mutations
total_gost<- gost(list("Genes with mutations" = rownames(maf_data_df)[which(maf_data_df$MutatedSamples >=1)]), 
                  organism = "mmusculus", 
                  ordered_query = FALSE, 
                  multi_query = FALSE, 
                  significant = TRUE, 
                  exclude_iea = FALSE,
                  measure_underrepresentation = FALSE, 
                  evcodes = TRUE,
                  user_threshold = 0.05, 
                  correction_method = "fdr",
                  domain_scope = "custom", 
                  custom_bg =custom_bg,
                  numeric_ns = "", 
                  sources = NULL, 
                  as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

ggplot(results[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size = 12, color ="black"),
        axis.text.y = element_text(size = 12, color ="black")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")


#############| ENRICHMENT
#|      * ALL mutations in the subset
subset <- c("Sample_9","Sample_5","Sample_4","Sample_14","Sample_8","Sample_2","Sample_7","Sample_17","Sample_25")
filtered_mutations_df_subset <- filtered_mutations_df[which(filtered_mutations_df$Samples %in% subset),]
maf_df <- filtered_mutations_df_subset %>%
  mutate(
    Tumor_Sample_Barcode = Samples,        # Map sample names
    Hugo_Symbol = SYMBOL,                  # Map gene names
    Chromosome = "Unknown",                # Placeholder if chromosome not available
    Start_Position = POS,                  # Use POS as start
    End_Position = POS,                    # Use POS as end for single-base mutations
    VAF = AF,
    Variant_Classification = Consequence_maftools,
    Variant_Type = "SNV",                  # Placeholder for mutation type
    Reference_Allele = "-",                # Placeholder if data not available
    Tumor_Seq_Allele2 = "-",               # Placeholder if data not available
    Protein_Change = NA_character_         # Placeholder if not available
  ) %>%
  dplyr::select(
    Tumor_Sample_Barcode, 
    Hugo_Symbol, Chromosome, Start_Position, End_Position,VAF,
    Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele2, Protein_Change
  )

#| DEFINING COLORS
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Set2')
names(vc_cols) = c(
  'Missense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Frame_Shift_Del',
  'Splice_Site',
  'In_Frame_Del',
  'Nonsense_Mutation'
)

maf_data <- read.maf(maf = maf_df)
maf_data_df <- as.data.frame(maf_data@gene.summary)
rownames(maf_data_df) <- maf_data_df$Hugo_Symbol
maf_data_df <- maf_data_df[,-c(1)]

#| INCLUDING ALLELE FREQUENCY
genes <- rownames(maf_data_df)[which(maf_data_df$MutatedSamples >=2)]
maf_data@gene.summary$Hugo_Symbol
aml_genes <- maf_data@gene.summary$Hugo_Symbol[c(1:length(genes))]
aml_genes_vaf = subsetMaf(maf = maf_data, 
                          genes = aml_genes, 
                          fields = "VAF", mafObj = FALSE)[,mean(VAF, na.rm = TRUE), Hugo_Symbol]

colnames(aml_genes_vaf)[2] = "VAF"
rownames(maf_data_df)[which(maf_data_df$AlteredSamples >=1)]

total_gost<- gost(list("Genes with mutations" = rownames(maf_data_df)[which(maf_data_df$MutatedSamples >=1)]), 
                  organism = "mmusculus", 
                  ordered_query = FALSE, 
                  multi_query = FALSE, 
                  significant = TRUE, 
                  exclude_iea = FALSE,
                  measure_underrepresentation = FALSE, 
                  evcodes = TRUE,
                  user_threshold = 0.05, 
                  correction_method = "fdr",
                  domain_scope = "custom", 
                  custom_bg =custom_bg,
                  numeric_ns = "", 
                  sources = NULL, 
                  as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")
ggplot(results[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size = 12, color ="black"),
        axis.text.y = element_text(size = 12, color ="black")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
################################################################################