################################################################################
###   CODEX2   ###
################################################################################

#| Source: https://htmlpreview.github.io/?https://github.com/yuchaojiang/CODEX2/blob/master/demo/CODEX2.html

#| In the case of mm39, I needed to change a bit the script using the source from:
#|    * https://github.com/yuchaojiang/CODEX2/tree/master/mouse

#| In the case of genome mm39, mm10 or mm9, they have not include in their pipeline the
#| mappability. I was trying to find it or compute it but the information in the UCSC
#| is not complete. I need to compute it myself using the code they suggest in here:
#|    * https://groups.google.com/g/codex2/c/x8b390WZNI0/m/-2UTtvLBAwAJ
#| However, CODEX2 does not require mappability directly as an input.  It is not strictly 
#| necessary to compute mappability yourself because the tool does not rely on mappability 
#| as a direct input. 
#| Mappability can be relevant when:
#|    * You are working with whole-genome sequencing (WGS) rather than targeted panels.
#|    * The regions of interest have a high degree of repetitiveness or low uniqueness (e.g., centromeres or segmental duplications).
#|    *You are designing or evaluating a new target panel.
#| In such cases, mappability information can help filter out low-quality regions where CNV detection might be unreliable.
#################################################################################


################################################################################
#| Functions
################################################################################
getgc.mm10 =function (ref) {
  seqlevelsStyle(ref)='UCSC'
  gc=rep(NA,length(ref))
  for(chr in unique(seqnames(ref))){
    message("Getting GC content for ", chr, sep = "")
    chr.index=which(as.matrix(seqnames(ref))==chr)
    ref.chr=IRanges(start= start(ref)[chr.index] , end = end(ref)[chr.index])
    if (chr == "X" | chr == "x" | chr == "chrX" | chr == "chrx") {
      chrtemp <- 20
    } else if (chr == "Y" | chr == "y" | chr == "chrY" | chr == "chry") {
      chrtemp <- 21
    } else if (chr == "M" | chr == "m" | chr == "chrM" | chr == {"chrm"}){
      chrtemp <- 22
    } else {
      chrtemp <- as.numeric(mapSeqlevels(as.character(chr), "NCBI")[1])
    }
    if (length(chrtemp) == 0) message("Chromosome cannot be found in NCBI database!")
    chrm <- unmasked(Mmusculus[[chrtemp]])
    seqs <- Views(chrm, ref.chr)
    af <- alphabetFrequency(seqs, baseOnly = TRUE, as.prob = TRUE)
    gc[chr.index] <- round((af[, "G"] + af[, "C"]) * 100, 2)
  }
  gc
}

#| This I modified myself from the human version
getgc_modifed <- function (ref, genome = NULL) {
  if(is.null(genome)){genome=BSgenome.Hsapiens.UCSC.hg19}
  gc=rep(NA,length(ref))
  for(chr in unique(seqnames(ref))){
    message("Getting GC content for chr ", chr, sep = "")
    chr.index=which(as.matrix(seqnames(ref))==chr)
    ref.chr=IRanges(start= start(ref)[chr.index] , end = end(ref)[chr.index])
    if (chr == "X" | chr == "x" | chr == "chrX" | chr == "chrx") {
      chrtemp <- 'X'
    } else if (chr == "Y" | chr == "y" | chr == "chrY" | chr == 
               "chry") {
      chrtemp <- 'Y'
    } else {
      chrtemp <- as.numeric(mapSeqlevels(as.character(chr), "NCBI")[1])
    }
    if (length(chrtemp) == 0) message("Chromosome cannot be found in NCBI Homo sapiens database!")
    chrm <- unmasked(genome[[paste('chr',chrtemp,sep='')]])
    seqs <- Views(chrm, ref.chr)
    af <- alphabetFrequency(seqs, baseOnly = TRUE, as.prob = TRUE)
    gc[chr.index] <- round((af[, "G"] + af[, "C"]) * 100, 2)
  }
  gc
}

#| Compute mappability by chromosome
computemapp<-function(ref.chr, L, genome, chr){
  genome.chr<-genome[[chr]]
  mapp.chr=rep(1,length(ref.chr))
  for(mappi in 1:length(mapp.chr)){
    cat(mappi,'\t')
    if(width(ref.chr)[mappi]>=L){
      dict = Views(genome.chr, start=seq((start(ref.chr))[mappi],(end(ref.chr))[mappi]-L+1,1), width=L)
    } else{
      dict = Views(genome.chr, start=start(ref.chr)[mappi]+round(width(ref.chr)[mappi]/2)-round(L/2), width=L)
    }
    if(sum(alphabetFrequency(dict, baseOnly=T)[,'other'])>0){
      mapp.chr[mappi]=0
    } else{
      pd = PDict(dict)
      ci=rep(0,length(pd))
      for(t in 1:length(genome)){
        res=matchPDict(pd,genome[[t]]); ci=ci+elementNROWS(res)
      }
      mapp.chr[mappi]=1/mean(ci)
    }
  }
  return(mapp.chr)
}

################################################################################


################################################################################
#| LIBRARY AND DATA
################################################################################
#| options(timeout = 600)  # Increase timeout to 10 minutes
#| BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggplot2))
suppressMessages(library(GenomicRanges))
suppressMessages(library(gprofiler2))
suppressMessages(library(IRanges))
suppressMessages(library(pheatmap))
suppressMessages(library(GenomicRanges))
suppressMessages(library(biomaRt))
suppressMessages(library(gprofiler2))
suppressMessages(library(viridis))
suppressMessages(library(GenVisR))
suppressMessages(library(absCNseq))
suppressMessages(library(dplyr))
suppressMessages(library(maftools))
suppressMessages(library(CODEX2))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(tidyr))
suppressMessages(library(dotenv))

#| For plots
theme_set(theme_classic())

#| Setting working dir
setwd("../")

# Configuration env file
load_dot_env(file = "../00_conf_env.env") 

#| Dataset from ensembl
genome <- "mmusculus_gene_ensembl"
version <- 102

#| Directory Results
dir.results <- paste0(Sys.getenv("R_BASE_DIR"),Sys.getenv("CNAs_CODEX2_DIR"), "Results/")

dir.create(paste0(dir.results,"ViolinPlots"))
dir.create(paste0(dir.results,"CODEX2_Metrics"))
dir.create(paste0(dir.results,"Heatmaps"))
dir.create(paste0(dir.results,"Histograms"))

#| Mappability scores
dir.mappability <- paste0(Sys.getenv("R_DATA_SHARED"),Sys.getenv("MAPPABILITY_DIR"))

#| Reading mouse information (containing gene type info)
ensembl <- useEnsembl(biomart = "genes", dataset = genome, version = version)
annotations <- getBM(attributes = c("gene_biotype", "external_gene_name", "chromosome_name", "start_position", "end_position"), mart = ensembl)
names(annotations) <- c("gene_type","gene_name", "CONTIG", "START", "END" )
annotations$CONTIG <- paste("chr",annotations$CONTIG, sep ="")

##########| Pre-processing

#| Extracting normal and tumor sample ID
sampname <- c(strsplit(Sys.getenv("TUMOR_PON_SAMPLES"), ",")[[1]] ,strsplit(Sys.getenv("NORMAL_PON_SAMPLES"), ",")[[1]]  )

#| Extracting normal sample ID
sampname_norm <-  strsplit(Sys.getenv("NORMAL_PON_SAMPLES"), ",")[[1]] 

#| Extracting tumor sample ID bam files
samples_tumor <- paste(strsplit(Sys.getenv("TUMOR_PON_SAMPLES"), ",")[[1]], "_mapped.rg.sorted.dup.bqsr.bam",sep ="")

#| Indicating CNAs directory
dirPath <- paste0(Sys.getenv("R_BASE_DIR"),Sys.getenv("MAPPED_DIRECTORY"))

#| BED file
bedFile <- paste0(Sys.getenv("R_DATA_SHARED"),Sys.getenv("BED_FILE"))

#| Dir with bam files
bamdir <- paste(dirPath,paste(sampname, "_mapped.rg.sorted.dup.bqsr.bam",sep =""), sep ="")

#| Creating bam object 
bambedObj <- getbambed(bamdir = bamdir, 
                       bedFile = bedFile, 
                       sampname = sampname, 
                       projectname = Sys.getenv("PROJECT_NAME"))

bamdir <- bambedObj$bamdir
sampname <- bambedObj$sampname
ref <- bambedObj$ref
projectname <- bambedObj$projectname
seqlevels(ref) <- paste0("chr", seqlevels(ref))

#| Getting mappability
mappability <- read.table(dir.mappability)

#| Saving mappabilty scores in a GRanger object
mappability_gr <- GRanges(
  seqnames = mappability$V1,                                       # Chromosome column
  ranges = IRanges(start = mappability$V2, end = mappability$V3),  # Start and end columns
  score = mappability$V5                                           # Mappability score
)
################################################################################




################################################################################
#| CODEX2 preprocessing and quality control
################################################################################

#| In cases when you do not know the mappability scores for the genome in study use:
#|    mapp <- DataFrame(data.frame(mapp = rep(1,dim(DataFrame(gc))[1]))) # Use in case of not finding the mappability scores for the genome analyzed
#|    values(ref) <- cbind(values(ref), DataFrame(gc,mapp))  


##########| Getting GC content and mappability scores
genome <- Mmusculus
gc <- getgc.mm10(ref)
values(ref) <- cbind(values(ref), DataFrame(gc))  

# Find overlaps between gr1 and gr2
overlaps <- findOverlaps(ref, mappability_gr)

# Extract matching `score` values from gr2
score_values <- rep(NA, length(ref))  # Default NA for non-overlapping regions
score_values[queryHits(overlaps)] <- mcols(mappability_gr)$score[subjectHits(overlaps)]

# Add the `score` column to gr1
mcols(ref)$mapp <- score_values
mcols(ref)$mapp[which(is.na(mcols(ref)$mapp))] <- 0

##########| Getting raw read depth
coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y
write.csv(Y, file = paste(projectname, '_coverage.csv', sep=''), quote = FALSE)

#| IF loaded from read.cvs do
Y <- read.csv(paste(projectname, '_coverage.csv', sep=''))
rownames(Y) <- Y$X
Y <- Y[,-c(1)]
Y <- as.matrix(Y)

#| visualizing Coverage Distribution
df_coverage <- data.frame(Y)
df_coverage$region <- rownames(df_coverage)
df_coverage <- gather(df_coverage, key = "Sample", value = "Value", -region)
df_coverage$Sample[which(df_coverage$Sample %in% sampname_norm)] <- gsub("Sample_","WT_",df_coverage$Sample[which(df_coverage$Sample %in% sampname_norm)])
df_coverage$Sample[which( !(df_coverage$Sample %in% sampname_norm))] <- gsub("Sample_","Tumour_",df_coverage$Sample[which(!(df_coverage$Sample %in% sampname_norm))])
df_coverage$Value_log2 <- log2(df_coverage$Value + 1)

ggplot(df_coverage, aes(x = reorder(Sample,-Value_log2), y = Value_log2, fill =Sample)) +
  geom_violin(trim=FALSE) +
  stat_summary(fun.data=mean_sdl, mult=1, 
               geom="pointrange", color="black") +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size =12, angle = 60, hjust=1, color ="black"),
        axis.text.y = element_text(size =12, color ="black"),
        legend.position = "none")+
  scale_fill_manual(values = rainbow(27))+
  ylab("Read depth")+
  xlab("")
ggsave(paste0(dir.results,"ViolinPlots/Read-depth_log2.pdf"), height = 4.5, width = 8)

ggplot(df_coverage, aes(x = reorder(Sample,-Value), y = Value, fill =Sample)) +
  geom_violin(trim=FALSE) +
  stat_summary(fun.data=mean_sdl, mult=1, 
               geom="pointrange", color="black") +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size =12, angle = 60, hjust=1, color ="black"),
        axis.text.y = element_text(size =12, color ="black"),
        legend.position = "none")+
  scale_fill_manual(values = rainbow(27))+
  ylab("Read depth")+
  xlab("")
ggsave(paste0(dir.results,"ViolinPlots/Read-depth.pdf"), height = 4.5, width = 8)


##########| Quality control 
qcObj <- qc(Y, 
            sampname, 
            ref, 
            cov_thresh = c(20, 3000),
            length_thresh = c(20, 2000), 
            mapp_thresh = 0.90,
            gc_thresh = c(20, 80))
Y_qc <- qcObj$Y_qc
sampname_qc <- qcObj$sampname_qc
ref_qc <- qcObj$ref_qc
qcmat <- qcObj$qcmat
gc_qc <- ref_qc$gc

##########| Running CODEX2
Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero,1,function(x){exp(1/length(x)*sum(log(x)))})
N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)
################################################################################



################################################################################
#| METHOD: CHROMOSOME-BY
################################################################################

finalcall.CBS_list <- list()
chromosomes <- unique(levels(seqnames(ref_qc)))
optimal_k_normalize <- c()

for (chr in chromosomes){
  
  print(chr)
  
  chr.index <- which(seqnames(ref_qc) %in% chr)
  sampname_norm_index <- which(colnames(Y_qc) %in% sampname_norm)

  normObj <- normalize_codex2_ns(Y_qc = Y_qc[chr.index,],
                                 gc_qc = gc_qc[chr.index], 
                                 K = 1:6, 
                                 norm_index = sampname_norm_index,
                                 N = N)
  Yhat.ns <- normObj$Yhat
  fGC.hat.ns <- normObj$fGC.hat
  beta.hat.ns <- normObj$beta.hat
  g.hat.ns <- normObj$g.hat 
  h.hat.ns <- normObj$h.hat
  AIC.ns <- normObj$AIC
  BIC.ns <- normObj$BIC
  RSS.ns <- normObj$RSS
  
  optimal_k_normalize <- c(optimal_k_normalize, which.max(BIC.ns))
  
  finalcall.CBS <- segmentCBS(Y_qc[chr.index,],  # recommended
                              Yhat.ns, 
                              optK = which.max(BIC.ns),
                              K = 1:6,
                              sampname_qc = colnames(Y_qc),
                              ref_qc = ranges(ref_qc)[chr.index],
                              chr = chr, 
                              lmax = 300, # If you want to call longer CNVs, you need to set the argument lmax to be larger. https://groups.google.com/g/codex2/c/r3Zc_zwTIWM/m/r3cJic9wAQAJ
                              mode = "fraction") #| for CNV detection in heterogeneous sample/tissue (e.g., somatic copy number changes in bulk cancer samples), use ‘fraction’ mode.
  
  finalcall.CBS_list[[chr]] <- finalcall.CBS
  
}

finalcall.CBS_list_df <- finalcall.CBS_list[[chromosomes[1]]]
for (i in 2:length(chromosomes)){
  finalcall.CBS_list_df <- rbind(finalcall.CBS_list_df, finalcall.CBS_list[[chromosomes[i]]])
}


write.table(finalcall.CBS_list_df, "finalcall.CBS_list_df_latent.txt", sep = "\t")



#| Chromosome 19 normalized coverage for Pten


gene_gr <- GRanges(
  seqnames = annotations$CONTIG,        # Chromosome/Contig information
  ranges = IRanges(start = annotations$START, end = annotations$END),  # Start and end positions
  gene_name = annotations$gene_name,
  gene_type = annotations$gene_type
)


# Find overlaps
overlaps <- findOverlaps(ref_qc, gene_gr)

overlap_data <- data.frame(
  chr = seqnames(ref_qc[queryHits(overlaps)]),
  index = queryHits(overlaps),
  start = start(ref_qc[queryHits(overlaps)]),
  end = end(ref_qc[queryHits(overlaps)]),
  gene_index = subjectHits(overlaps),
  gene_name = mcols(gene_gr[subjectHits(overlaps)])$gene_name,
  gene_type = mcols(gene_gr[subjectHits(overlaps)])$gene_type,
  gc_content = mcols(ref_qc[queryHits(overlaps)])$gc,
  mapp = mcols(ref_qc[queryHits(overlaps)])$mapp
)





overlap_data[which(overlap_data$gene_name == "Pten"),]

length(Y_qc[as.numeric(overlap_data$index[which(overlap_data$chr == "chr19")])])


finalcall.CBS_filtered_df <- finalcall.CBS_list[["chr19"]]
#| Create CNA GRanges object
cna_gr <- GRanges(
  seqnames = finalcall.CBS_filtered_df$chr,        
  ranges = IRanges(start = finalcall.CBS_filtered_df$st_bp, end =finalcall.CBS_filtered_df$ed_bp),
  copynumber = finalcall.CBS_filtered_df$copy_no,
  cnv = finalcall.CBS_filtered_df$cnv,
  length_kb =finalcall.CBS_filtered_df$length_kb,
  lratio = finalcall.CBS_filtered_df$lratio,
  sample_id = finalcall.CBS_filtered_df$sample_name,
  raw_cov = finalcall.CBS_filtered_df$raw_cov,
  norm_cov = finalcall.CBS_filtered_df$norm_cov,
  mBIC = finalcall.CBS_filtered_df$mBIC,
  st_exon = finalcall.CBS_filtered_df$st_exon,
  ed_exon = finalcall.CBS_filtered_df$ed_exon
)

#| Create Gene GRanges object
gene_gr <- GRanges(
  seqnames = annotations$CONTIG,        # Chromosome/Contig information
  ranges = IRanges(start = annotations$START, end = annotations$END),  # Start and end positions
  gene_name = annotations$gene_name,
  gene_type = annotations$gene_type
)


#| Find overlaps
overlaps <- findOverlaps(cna_gr, gene_gr)

overlap_data <- data.frame(
  chr = seqnames(cna_gr[queryHits(overlaps)]),
  cna_index = queryHits(overlaps),
  gene_index = subjectHits(overlaps),
  cna_width = width(cna_gr[queryHits(overlaps)]),
  gene_name = mcols(gene_gr[subjectHits(overlaps)])$gene_name,
  gene_type = mcols(gene_gr[subjectHits(overlaps)])$gene_type,
  cna_start = start(cna_gr[queryHits(overlaps)]),
  cna_end = end(cna_gr[queryHits(overlaps)]),
  gene_start = start(gene_gr[subjectHits(overlaps)]),
  gene_end = end(gene_gr[subjectHits(overlaps)]),
  copynumber = mcols(cna_gr[queryHits(overlaps)])$copynumber,
  type = mcols(cna_gr[queryHits(overlaps)])$cnv,
  lratio = mcols(cna_gr[queryHits(overlaps)])$lratio,
  sample_id = mcols(cna_gr[queryHits(overlaps)])$sample_id,
  length_kb = mcols(cna_gr[queryHits(overlaps)])$length_kb,
  raw_cov = mcols(cna_gr[queryHits(overlaps)])$raw_cov,
  norm_cov = mcols(cna_gr[queryHits(overlaps)])$norm_cov,
  mBIC = mcols(cna_gr[queryHits(overlaps)])$mBIC,
  st_exon = mcols(cna_gr[queryHits(overlaps)])$st_exon,
  ed_exon = mcols(cna_gr[queryHits(overlaps)])$ed_exon
)

overlap_data[which(overlap_data$gene_name == "Pten"),c("gene_name", "copynumber", "lratio", "raw_cov","sample_id","type")]
