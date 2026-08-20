################################################################################
###   CODEX2 DATA ANALYSIS OF AC-82_WESmouse  ###
################################################################################

#| Why sharing starting points or fragments across samples could be expected: 

#| Segmentation Algorithms:
#|    * CODEX2 uses segmentation to identify breakpoints in CNAs. If multiple samples have similar 
#|    coverage patterns at specific loci, the algorithm might segment these regions at the same positions, 
#|    resulting in shared starting or ending points.

#| Batch Effects:
#|    * CODEX2 relies on normalized coverage, which can be influenced by batch effects or shared biases 
#|    in the data (e.g., library prep or sequencing). These shared effects may lead to systematic similarities 
#|    in CNAs across samples.

#| Biological Relevance:
#|    * As with any CNA detection method, biologically recurrent CNAs (e.g., at oncogene loci or fragile sites) 
#|    will naturally show shared breakpoints across samples.

#| Mouse-Specific Features:
#|    * If you are working with inbred mice, the uniformity of their genomes may lead to CNAs with consistent start/end 
#|    points due to inherited structural variations.

#| Is This Expected in CODEX2?
#|    Yes, shared CNA starting points are expected under certain circumstances:
#|      - Recurrent CNAs: These are biologically relevant and occur consistently across samples due to their involvement in disease pathways.
#|      - Algorithm Behavior: CODEX2 may identify similar breakpoints across samples because of shared signal patterns in the coverage data.
#|    However, unexpected shared starting points can also arise:
#|      - Pipeline Artifacts: Segmentation algorithms sometimes favor specific regions due to low complexity or repetitive sequences.
#|      - Insufficient Normalization: If CODEX2's GC-content normalization or latent factor adjustment isn’t fully accounting for biases, artifacts can emerge.

#| Given that developers from CODEX2 haven't previously computed the mappability for mm39,
#| they recommend to run the program anyway by setting to 1 the vector of mappability,
#| in that way, the program does not get any error because of this step. The only
#| inconvenient is that we will not be able to filter by low mappability regions. see https://groups.google.com/g/codex2/c/TZOXXdYHnwo/m/JEX6cOejAAAJ
#| However, in cases where you have high sequencing depth, even if mappability is low 
#| in certain regions, you might still obtain enough reads from other regions to ensure 
#| accurate results. The higher the coverage, the more likely you are to get a sufficient 
#| number of unique reads, even from repetitive areas.
#| For targeted sequencing approaches, where you're sequencing a predefined set of regions 
#| (such as the exome or specific gene panels), mappability may be less of an issue because 
#| the regions of interest are usually well-characterized and designed to avoid problematic, 
#| repetitive sequences. The focus here is on accurately capturing the target regions, and 
#| mappability might be less critical as long as those regions are covered.


#| For cancer research projects related to CODEX2, the have applied the filter of copy number < 1.7
#| for deletions and > 2.3 for duplications https://groups.google.com/g/codex2/c/1vTscWGDg18/m/rMV8LayDAwAJ

#| Justification of the enrichment: If you want to retain enrichment results for pathways with 
#| few genes, using the whole genome as the background is reasonable. This maximizes the number 
#| of possible annotations and increases sensitivity for small gene sets.
################################################################################+


################################################################################
#| LIBRARY AND DATA
################################################################################
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
suppressMessages(library(tidyr))
suppressMessages(library(colorspace))
suppressMessages(library(ggrepel))
suppressMessages(library(patchwork))
suppressMessages(library(dotenv))

#| For plots
theme_set(theme_classic())

#| Setting working dir
setwd("../")

# Configuration env file
load_dot_env(file = "../00_conf_env.env") 

#| Dataset from ensembl
genome <- Sys.getenv("DATASET_ENSEMBL")
version <- as.numeric(Sys.getenv("VERSION_ENSEMBL"))

#| Directory Results
dir.results <- paste0(Sys.getenv("R_BASE_DIR"),Sys.getenv("CNAs_CODEX2_DIR"), "Results/")

#| Extracting normal and tumor sample ID
sampname <- c(strsplit(Sys.getenv("TUMOR_PON_SAMPLES"), ",")[[1]] ,strsplit(Sys.getenv("NORMAL_PON_SAMPLES"), ",")[[1]]  )

#| Extracting normal sample ID
sampname_norm <-  strsplit(Sys.getenv("NORMAL_PON_SAMPLES"), ",")[[1]] 

#| Extracting tumor sample ID bam files
samples_tumor <- paste(strsplit(Sys.getenv("TUMOR_PON_SAMPLES"), ",")[[1]], "_mapped.rg.sorted.dup.bqsr.bam",sep ="")

#| Reading mouse information (containing gene type info)
ensembl <- useEnsembl(biomart = "genes", dataset = genome, version = version)
annotations <- getBM(attributes = c("gene_biotype", "external_gene_name", "chromosome_name", "start_position", "end_position"), mart = ensembl)
names(annotations) <- c("gene_type","gene_name", "CONTIG", "START", "END" )
annotations$CONTIG <- paste0("chr",annotations$CONTIG)

#| Custom background for enrichment analysis
custom_bg <- unique(annotations$gene_name[which(annotations$gene_type == "protein_coding")])

#| Reading the CNAs from CODEX2 results
finalcall.CBS_list_df <- read.table(paste0(Sys.getenv("R_BASE_DIR"),Sys.getenv("CNAs_CODEX2_DIR"), "finalcall.CBS_list_df_latent.txt"), sep ="\t")

#| Lkb distribution
ggplot(finalcall.CBS_list_df, aes(x = length_kb ))+
  geom_histogram(color="black",binwidth=200)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.y = element_text(size = 12,color ="black"),
        axis.text.x = element_text(size = 12,color ="black"),
        plot.title=element_text(size=12)) +
  ylab("Counts")+
  xlab("CNA Length")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_length_kb.pdf"), height = 4.5, width = 5.5)

#| lratio distribution
ggplot(finalcall.CBS_list_df, aes(x = lratio))+
  geom_histogram(color="black",binwidth=20)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.y = element_text(size = 12,color ="black"),
        axis.text.x = element_text(size = 12,color ="black"),
        plot.title=element_text(size=12)) +
  ylab("Counts")+
  xlab("lratio")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_lratio.pdf"), height = 4.5, width = 5.5)

#| mBIC distribution
ggplot(finalcall.CBS_list_df, aes(x = mBIC))+
  geom_histogram(color="black",binwidth=200)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.y = element_text(size = 12,color ="black"),
        axis.text.x = element_text(size = 12,color ="black"),
        plot.title=element_text(size=12)) +
  ylab("Counts")+
  xlab("mBIC")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_mBIC.pdf"), height = 4.5, width = 5.5)

#| Copy number distribution
ggplot(finalcall.CBS_list_df, aes(x = copy_no, fill =cnv))+
  geom_histogram(binwidth=0.05, color = "black")+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.y = element_text(size = 12,color ="black"),
        axis.text.x = element_text(size = 12,color ="black"),
        plot.title=element_text(size=12)) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("Copy number")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_copynumber.pdf"), height = 4.5, width = 5.5)

#| Raw coverage distribution
ggplot(finalcall.CBS_list_df, aes(x = raw_cov))+
  geom_histogram(binwidth=500, color = "black")+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("Raw coverage")+
  ylim(0,510)
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_Raw-Coverage.pdf"), height = 4.5, width = 5.5)

#| Mean Raw coverage distribution
ggplot(finalcall.CBS_list_df, aes(x = norm_cov))+
  geom_histogram(binwidth=500, color = "black")+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("Normalized Raw Coverage")+
  ylim(0,510)
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_Normalized-Raw-Coverage.pdf"), height = 4.5, width = 5.5)

#| Raw vs normalized coverage
cor <- cor.test(finalcall.CBS_list_df$norm_cov, finalcall.CBS_list_df$raw_cov, method ="pearson", exact = FALSE)
p <- round(cor$p.value,abs(floor(log10(abs(cor$p.value)))))
r <- round(as.numeric(cor$estimate[[1]]),2)
ggplot(finalcall.CBS_list_df, aes(x = norm_cov, y =raw_cov))+
  geom_point(color="black",size=2, alpha =0.2)+
  geom_smooth(method = "lm", color ="green2") +
  theme(text=element_text(size=14,  family="sans"),
        axis.text.y = element_text(size = 12,color ="black"),
        axis.text.x = element_text(size = 12,color ="black"),
        plot.title=element_text(size=12)) +
  ylab("Raw coverage")+
  xlab("Normalized Raw Coverage")+
  ggtitle(paste("r = ",r,", p < ", format(p, digits =2, scientific = T), sep=""))
ggsave(paste0(dir.results,"CODEX2_Metrics/Scatter_Normalized-Raw-Coverage.pdf"), height = 4.5, width = 5)

#| Raw vs normalized coverage
cor <- cor.test(finalcall.CBS_list_df$mBIC, finalcall.CBS_list_df$raw_cov, method ="pearson", exact = FALSE)
p <- round(cor$p.value,abs(floor(log10(abs(cor$p.value)))))
r <- round(as.numeric(cor$estimate[[1]]),2)
ggplot(finalcall.CBS_list_df, aes(x = mBIC, y =raw_cov))+
  geom_point(color="black",size=2, alpha =0.2)+
  geom_smooth(method = "lm", color ="green2") +
  theme(text=element_text(size=14,  family="sans"),
        axis.text.y = element_text(size = 12,color ="black"),
        axis.text.x = element_text(size = 12,color ="black"),
        plot.title=element_text(size=12)) +
  ylab("Raw coverage")+
  xlab("mBIC")+
  ggtitle(paste("r = ",r,", p < ", format(p, digits =2, scientific = T), sep=""))
ggsave(paste0(dir.results,"CODEX2_Metrics/Scatter_mBIC-Raw-Coverage.pdf"), height = 4.5, width = 5)


ggplot(finalcall.CBS_list_df, aes(x =copy_no , y =mBIC ))+
  geom_point(color="black",size=2, alpha =0.2)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.y = element_text(size = 12,color ="black"),
        axis.text.x = element_text(size = 12,color ="black"),
        plot.title=element_text(size=12)) +
  geom_hline(yintercept = 5000, color ="red")+
  geom_vline(xintercept = 1.7, color ="green")+
  geom_vline(xintercept = 2.3, color ="green")+
  xlab("copy_no")+
  ylab("mBIC")
ggsave(paste0(dir.results,"CODEX2_Metrics/Scatter_mBIC-copy_no.pdf"), height = 4.5, width = 5)

ggplot(finalcall.CBS_list_df, aes(x =copy_no , y =length_kb ))+
  geom_point(color="black",size=2, alpha =0.2)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.y = element_text(size = 12,color ="black"),
        axis.text.x = element_text(size = 12,color ="black"),
        plot.title=element_text(size=12)) +
  geom_hline(yintercept = 200, color ="red")+
  geom_vline(xintercept = 1.7, color ="green")+
  geom_vline(xintercept = 2.3, color ="green")+
  xlab("copy_no")+
  ylab("Length")
ggsave(paste0(dir.results,"CODEX2_Metrics/Scatter_length-copy_no.pdf"), height = 4.5, width = 5)

ggplot(finalcall.CBS_list_df, aes(x =copy_no , y =norm_cov ))+
  geom_point(color="black",size=2, alpha =0.2)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.y = element_text(size = 12,color ="black"),
        axis.text.x = element_text(size = 12,color ="black"),
        plot.title=element_text(size=12)) +
  geom_vline(xintercept = 1.7, color ="green")+
  geom_vline(xintercept = 2.3, color ="green")+
  xlab("copy_no")+
  ylab("Normalized Raw Coverage")
ggsave(paste0(dir.results,"CODEX2_Metrics/Scatter_Normalized-Raw-Coverage-copy_no.pdf"), height = 4.5, width = 5)

ggplot(finalcall.CBS_list_df, aes(x =mBIC , y =lratio ))+
  geom_point(color="black",size=2, alpha =0.2)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.y = element_text(size = 12,color ="black"),
        axis.text.x = element_text(size = 12,color ="black"),
        plot.title=element_text(size=12)) +
  xlab("copy_no")+
  ylab("Normalized Raw Coverage")
#ggsave(paste0(dir.results,"CODEX2_Metrics/Scatter_Normalized-Raw-Coverage-copy_no.pdf"), height = 4.5, width = 5)




#| How is Pten when no filters are applied?

sapply(finalcall.CBS_list_df, length)

# Create CNA GRanges object
cna_gr <- GRanges(
  seqnames = finalcall.CBS_list_df$chr,        
  ranges = IRanges(start = finalcall.CBS_list_df$st_bp, end =finalcall.CBS_list_df$ed_bp),
  copynumber = finalcall.CBS_list_df$copy_no,
  cnv = finalcall.CBS_list_df$copy_no,
  length_kb =finalcall.CBS_list_df$length_kb,
  lratio = finalcall.CBS_list_df$lratio,
  sample_id = finalcall.CBS_list_df$sample_name,
  raw_cov = finalcall.CBS_list_df$raw_cov,
  norm_cov = finalcall.CBS_list_df$norm_cov,
  mBIC = finalcall.CBS_list_df$mBIC,
  st_exon = finalcall.CBS_list_df$st_exon,
  ed_exon = finalcall.CBS_list_df$ed_exon
)

# Create Gene GRanges object
gene_gr <- GRanges(
  seqnames = annotations$CONTIG,        # Chromosome/Contig information
  ranges = IRanges(start = annotations$START, end = annotations$END),  # Start and end positions
  gene_name = annotations$gene_name,
  gene_type = annotations$gene_type
)


# Find overlaps
overlaps <- findOverlaps(cna_gr, gene_gr)

overlap_data_no_filters <- data.frame(
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
  copynumber = mcols(cna_gr[queryHits(overlaps)])$cnv,
  lratio = mcols(cna_gr[queryHits(overlaps)])$lratio,
  sample_id = mcols(cna_gr[queryHits(overlaps)])$sample_id,
  length_kb = mcols(cna_gr[queryHits(overlaps)])$length_kb ,
  raw_cov = mcols(cna_gr[queryHits(overlaps)])$raw_cov,
  norm_cov = mcols(cna_gr[queryHits(overlaps)])$norm_cov,
  mBIC = mcols(cna_gr[queryHits(overlaps)])$mBIC,
  st_exon = mcols(cna_gr[queryHits(overlaps)])$st_exon,
  ed_exon = mcols(cna_gr[queryHits(overlaps)])$ed_exon
)

#| Filtering places with no gene associated
overlap_data_no_filters <- overlap_data_no_filters[which(overlap_data_no_filters$gene_name != ""),]

#| Filtering genes appearing in the panel of normal
overlap_data_no_filters <- overlap_data_no_filters[which( !(overlap_data_no_filters$gene_name %in%  overlap_data_no_filters$gene_name[which(overlap_data_no_filters$sample_id %in% sampname_norm)]) ),]

overlap_data_no_filters_df <- overlap_data_no_filters %>%
  group_by(gene_name) %>%
  mutate(count_genes = n()) %>%
  ungroup()

overlap_data_no_filters_df <- overlap_data_no_filters_df %>%
  group_by(sample_id) %>%
  mutate(count_samples = n()) %>%
  ungroup()

overlap_data_no_filters_df[which(overlap_data_no_filters_df$gene_name == "Pten"),c("gene_type","sample_id","cna_width","copynumber","raw_cov","mBIC", "lratio")]
################################################################################


################################################################################
###| CODEX2 metrics evaluation 1
################################################################################

Lkb_cna <- 300
Exons <- 2
mBIC <- 4000
lratio <- 35
raw_cov <- 5
copy_no_del <- 1.6
copy_no_dup <- 2.4

#| Maximun size of the CNA & Ratio between the CNA and the exon
finalcall.CBS_filtered_df <- finalcall.CBS_list_df[which( (finalcall.CBS_list_df$length_kb<=Lkb_cna) & (finalcall.CBS_list_df$length_kb/(finalcall.CBS_list_df$ed_exon-finalcall.CBS_list_df$st_exon+1)<50) ),]

#| Number of exones in the CNA regions |  Minimum lratio
finalcall.CBS_filtered_df <- finalcall.CBS_filtered_df[which(  (finalcall.CBS_filtered_df$lratio>lratio) |  ((finalcall.CBS_filtered_df$ed_exon-finalcall.CBS_filtered_df$st_exon)>=Exons)),]

#| mBIC threshold
finalcall.CBS_filtered_df <- finalcall.CBS_filtered_df[which(finalcall.CBS_filtered_df$mBIC < mBIC),]

#| Copy-number filtering
finalcall.CBS_filtered_df <-finalcall.CBS_filtered_df[which( (finalcall.CBS_filtered_df$copy_no > copy_no_dup) |  (finalcall.CBS_filtered_df$copy_no < copy_no_del) ),]

#| Raw coverage
finalcall.CBS_filtered_df <- finalcall.CBS_filtered_df[which(finalcall.CBS_filtered_df$raw_cov > raw_cov),]

#| Raw vs normalized coverage
ggplot(finalcall.CBS_filtered_df, aes(x = norm_cov, y =raw_cov))+
  geom_point(color="black",size=2, alpha =0.2)+
  geom_smooth(method = "lm", color ="green") +
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  ylab("raw_cov")+
  xlab("norm_cov")
ggsave(paste0(dir.results,"CODEX2_Metrics/Scatter_Normalized-Raw-Coverage.pdf"), height = 4.5, width = 5.5)

ggplot(finalcall.CBS_filtered_df, aes(x = mBIC, y =norm_cov))+
  geom_point(color="black",size=1)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  ylab("norm_cov")+
  xlab("mBIC")

ggplot(finalcall.CBS_filtered_df, aes(x = mBIC, y =raw_cov))+
  geom_point(color="black",size=1)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  ylab("raw_cov")+
  xlab("mBIC")

length(finalcall.CBS_filtered_df$length_kb[which(finalcall.CBS_filtered_df$length_kb < 1)])
length(finalcall.CBS_filtered_df$length_kb[which(finalcall.CBS_filtered_df$length_kb > 1)])

#| Lkb distribution
ggplot(finalcall.CBS_filtered_df, aes(x = length_kb))+
  geom_histogram(color="black",binwidth=0.1)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("length_kb")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_length_filtered_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

#| lratio distribution
ggplot(finalcall.CBS_filtered_df, aes(x = lratio))+
  geom_histogram(color="black",binwidth=5)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("lratio")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_lratio_filtered_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

#| mBIC distribution
ggplot(finalcall.CBS_filtered_df, aes(x = mBIC))+
  geom_histogram(color="black",binwidth=10)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("mBIC")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_mBIC_filtered_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

#| Raw counts distribution
ggplot(finalcall.CBS_filtered_df, aes(x = raw_cov))+
  geom_histogram(color="black",binwidth=10)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("raw_cov ")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_Raw-Coverage_filtered_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

#| Normalized raw counts distribution
ggplot(finalcall.CBS_filtered_df, aes(x = norm_cov))+
  geom_histogram(color="black",binwidth=10)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("norm_cov")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_Normalized-Raw-Coverage_filtered_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)




#| Create CNA GRanges object
cna_gr <- GRanges(
  seqnames = finalcall.CBS_filtered_df$chr,        
  ranges = IRanges(start = finalcall.CBS_filtered_df$st_bp, end =finalcall.CBS_filtered_df$ed_bp),
  copynumber = finalcall.CBS_filtered_df$copy_no,
  cnv = finalcall.CBS_filtered_df$cnv,
  length_kb  =finalcall.CBS_filtered_df$length_kb ,
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
  length_kb  = mcols(cna_gr[queryHits(overlaps)])$length_kb ,
  raw_cov = mcols(cna_gr[queryHits(overlaps)])$raw_cov,
  norm_cov = mcols(cna_gr[queryHits(overlaps)])$norm_cov,
  mBIC = mcols(cna_gr[queryHits(overlaps)])$mBIC,
  st_exon = mcols(cna_gr[queryHits(overlaps)])$st_exon,
  ed_exon = mcols(cna_gr[queryHits(overlaps)])$ed_exon
)

#| Filtering places with no gene associated
overlap_data <- overlap_data[which(overlap_data$gene_name != ""),]

#| Filtering genes appearing in the panel of normal
overlap_data <- overlap_data[which( !(overlap_data$gene_name %in%  overlap_data$gene_name[which(overlap_data$sample_id %in% sampname_norm)]) ),]

overlap_data_df <- overlap_data %>%
  group_by(gene_name) %>%
  mutate(count_genes = n()) %>%
  ungroup()

overlap_data_df <- overlap_data_df %>%
  group_by(sample_id,gene_name) %>%
  mutate(count_genes_sample = n()) %>%
  ungroup()

overlap_data_df <- overlap_data_df %>%
  group_by(sample_id) %>%
  mutate(count_samples = n()) %>%
  ungroup()

overlap_data_df <- overlap_data_df %>%
  group_by(gene_name) %>%
  mutate(count_genes_sample_all = sum(count_genes_sample)) %>%
  ungroup()

writexl::write_xlsx(overlap_data_df, paste0(dir.results,"CNVs_results_filtered_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,".xlsx", sep = ""))

gene <- "Pten"
overlap_data_df$type[which(overlap_data_df$gene_name == gene)]
overlap_data_df$sample_id[which(overlap_data_df$gene_name == gene)]

#########| Plotting Histogram
ggplot(overlap_data, aes(x = copynumber, fill = type))+
  geom_histogram(color="black",binwidth=0.1)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("Copy number")
ggsave(paste0(dir.results,"Histograms/Copy_Number_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)


#########| Plotting barplot
overlap_data_dup_del <- overlap_data %>%
  group_by(sample_id,type) %>%
  summarise(count = n())

ggplot(overlap_data_dup_del, aes(x = reorder(sample_id,-count), y = count, fill = type ))+
  geom_bar(stat = "identity")+
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(color ="black")) +
  xlab("Samples")+
  ylab("Number of genes")+
  scale_fill_manual(values =c( "blue","red"))
ggsave(paste0(dir.results,"Histograms/Dup_Del_Sample_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 6)


############| Chr plot

overlap_data_df$mid_point <- (overlap_data_df$cna_start + overlap_data_df$cna_end)/2

ggplot(overlap_data_df, aes(x = mid_point, y = copynumber , color = copynumber )) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  scale_color_gradient2(low = "blue", 
                        high = "red", 
                        mid = "gray",
                        midpoint = 2) +
  geom_hline(yintercept = 0, alpha =0.5, color ="gray" ) +
  geom_vline(xintercept = seq(1, 100, by = 10), 
             color = "black", linetype = "dotted", alpha = 0.2)+
  geom_hline(yintercept = 2, color ="black") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    strip.text = element_text(size = 8),
    legend.position = "none"
  )  +
  labs(y = "Copy Number", x = "Genomic Position")
ggsave(paste0(dir.results,"Heatmaps/CNAs_Chromosome_Position-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 12)


#| Plotting CNVs by chromosomes
unique_cna <- overlap_data_df %>% distinct(cna_index, .keep_all = TRUE)

new_df <- unique_cna[,c("chr","cna_start", "cna_end", "copynumber", "sample_id")]
names(new_df) <- c("chromosome", "start", "end", "segmean", "sample")
new_df$chromosome <- gsub("chr", "", new_df$chromosome)
new_df$segmean <- as.integer(new_df$segmean)

cnFreq(new_df, genome="mm10",
       CN_Loss_colour="blue",
       CN_Gain_colour ="red",
       CN_low_cutoff = 1.7,
       CN_high_cutoff = 2.3) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        text=element_text(size=9,  family="sans"),
        axis.title.y = element_text(face = "bold", size =12),
        axis.title.x = element_text(size =12))
ggsave(paste0(dir.results,"Heatmaps/cnFreq_Proportion_CopyNumber_Gain_Loss-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 15)

cnSpec(
  x = new_df, 
  genome = "mm10", 
  plot_title = "Global Copy Number Alterations (CNA)",
  CN_Loss_colour = "blue", 
  CN_Gain_colour = "red"
)+theme_minimal() +  # Apply a minimal theme
  theme(
    legend.title = element_blank(),  # Remove legend title
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    strip.text.y=element_text(angle=0)
  )
ggsave(paste0(dir.results,"Heatmaps/cnSpec_Gain_Loss-Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 5.5, width = 9)


#########|  HEATMAP top candidates More than 6
overlap_data_df_f <- overlap_data_df[which(overlap_data_df$gene_type == "protein_coding"),]
overlap_data_df_f <- overlap_data_df_f[which( (overlap_data_df_f$count_genes > 5) & (overlap_data_df_f$count_genes_sample_all == overlap_data_df_f$count_genes ) ),]
overlap_data_df_f <- overlap_data_df_f %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon,count_genes,count_genes_sample_all) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )
heatmap_plot <- ggplot(overlap_data_df_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
overlap_data_df_f$freq_counts <- (overlap_data_df_f$count_genes/21)*100
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- paste(overlap_data_df_f$freq_counts_2, "%")
overlap_data_df_f$freq_counts <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)

#| Unique
unique_gene <- overlap_data_df_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/TOP_Candidates-6Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 5, width = 8)


#########|  HEATMAP top candidates More than 5
overlap_data_df_f <- overlap_data_df[which(overlap_data_df$gene_type == "protein_coding"),]
overlap_data_df_f <- overlap_data_df_f[which( (overlap_data_df_f$count_genes > 4) & (overlap_data_df_f$count_genes_sample_all == overlap_data_df_f$count_genes ) ),]
overlap_data_df_f <- overlap_data_df_f %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon,count_genes,count_genes_sample_all) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )
heatmap_plot <- ggplot(overlap_data_df_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
overlap_data_df_f$freq_counts <- (overlap_data_df_f$count_genes/21)*100
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- paste(overlap_data_df_f$freq_counts_2, "%")
overlap_data_df_f$freq_counts <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)

#| Unique
unique_gene <- overlap_data_df_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/TOP_Candidates-5Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 6.5, width = 8)


#########|  HEATMAP top candidates More than 4
overlap_data_df_f <- overlap_data_df[which(overlap_data_df$gene_type == "protein_coding"),]
overlap_data_df_f <- overlap_data_df_f[which( (overlap_data_df_f$count_genes > 3) & (overlap_data_df_f$count_genes_sample_all == overlap_data_df_f$count_genes ) ),]
overlap_data_df_f <- overlap_data_df_f %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon,count_genes,count_genes_sample_all) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )
heatmap_plot <- ggplot(overlap_data_df_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
overlap_data_df_f$freq_counts <- (overlap_data_df_f$count_genes/21)*100
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- paste(overlap_data_df_f$freq_counts_2, "%")
overlap_data_df_f$freq_counts <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)

#| Unique
unique_gene <- overlap_data_df_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/TOP_Candidates-4Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 6.5, width = 8)


#########|  HEATMAP top candidates More than 3
overlap_data_df_f <- overlap_data_df[which(overlap_data_df$gene_type == "protein_coding"),]
overlap_data_df_f <- overlap_data_df_f[which( (overlap_data_df_f$count_genes > 2) & (overlap_data_df_f$count_genes_sample_all == overlap_data_df_f$count_genes ) ),]
overlap_data_df_f <- overlap_data_df_f %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon,count_genes,count_genes_sample_all) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )
heatmap_plot <- ggplot(overlap_data_df_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
overlap_data_df_f$freq_counts <- (overlap_data_df_f$count_genes/21)*100
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- paste(overlap_data_df_f$freq_counts_2, "%")
overlap_data_df_f$freq_counts <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)

#| Unique
unique_gene <- overlap_data_df_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/TOP_Candidates-3Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 9, width = 10)


#########| HEATMAP WITH ONLY HUMAN HOMOLOG
homolog_info <- getBM(
  attributes = c("external_gene_name",  # Mouse Ensembl ID (for merging)
                 "hsapiens_homolog_associated_gene_name", 
                 "hsapiens_homolog_ensembl_gene", 
                 "hsapiens_homolog_orthology_type"),
  mart = ensembl
)
names(homolog_info) <- c(c("gene_name",  # Mouse Ensembl ID (for merging)
                           "hsapiens_homolog_associated_gene_name", 
                           "hsapiens_homolog_ensembl_gene", 
                           "hsapiens_homolog_orthology_type"))
annotations_homolog <- merge(annotations, homolog_info, by = "gene_name")
annotations_homolog <- annotations_homolog[which(annotations_homolog$hsapiens_homolog_associated_gene_name !=""),]

#| Selecting only homologs
overlap_data_homolog_df <- overlap_data_df[which(overlap_data_df$gene_name %in% unique(annotations_homolog$gene_name)),]


#########|  HEATMAP top candidates More than 4
overlap_data_homolog_df_f <- overlap_data_homolog_df[which(overlap_data_homolog_df$gene_type == "protein_coding"),]
overlap_data_homolog_df_f <- overlap_data_homolog_df_f[which( (overlap_data_homolog_df_f$count_genes > 3) & (overlap_data_homolog_df_f$count_genes_sample_all == overlap_data_homolog_df_f$count_genes ) ),]
overlap_data_homolog_df_f <- overlap_data_homolog_df_f %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon,count_genes,count_genes_sample_all) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )
heatmap_plot <- ggplot(overlap_data_homolog_df_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
overlap_data_homolog_df_f$freq_counts <- (overlap_data_homolog_df_f$count_genes/21)*100
overlap_data_homolog_df_f$freq_counts <- as.numeric(overlap_data_homolog_df_f$freq_counts)
overlap_data_homolog_df_f$freq_counts_2 <- round(overlap_data_homolog_df_f$freq_counts)
overlap_data_homolog_df_f$freq_counts_2 <- paste(overlap_data_homolog_df_f$freq_counts_2, "%")
overlap_data_homolog_df_f$freq_counts <- round(overlap_data_homolog_df_f$freq_counts)
overlap_data_homolog_df_f$freq_counts <- as.numeric(overlap_data_homolog_df_f$freq_counts)

#| Unique
unique_gene <- overlap_data_homolog_df_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/HOMOLOG-TOP_Candidates-4Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 6.5, width = 8)


#########|  HEATMAP top candidates More than 3
overlap_data_homolog_df_f <- overlap_data_homolog_df[which(overlap_data_homolog_df$gene_type == "protein_coding"),]
overlap_data_homolog_df_f <- overlap_data_homolog_df_f[which( (overlap_data_homolog_df_f$count_genes > 2) & (overlap_data_homolog_df_f$count_genes_sample_all == overlap_data_homolog_df_f$count_genes ) ),]
overlap_data_homolog_df_f <- overlap_data_homolog_df_f %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon,count_genes,count_genes_sample_all) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )
heatmap_plot <- ggplot(overlap_data_homolog_df_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
overlap_data_homolog_df_f$freq_counts <- (overlap_data_homolog_df_f$count_genes/21)*100
overlap_data_homolog_df_f$freq_counts <- as.numeric(overlap_data_homolog_df_f$freq_counts)
overlap_data_homolog_df_f$freq_counts_2 <- round(overlap_data_homolog_df_f$freq_counts)
overlap_data_homolog_df_f$freq_counts_2 <- paste(overlap_data_homolog_df_f$freq_counts_2, "%")
overlap_data_homolog_df_f$freq_counts <- round(overlap_data_homolog_df_f$freq_counts)
overlap_data_homolog_df_f$freq_counts <- as.numeric(overlap_data_homolog_df_f$freq_counts)

#| Unique
unique_gene <- overlap_data_homolog_df_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/HOMOLOG-TOP_Candidates-3Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 8.5, width = 9)


#########| HEATMAP WITH DATA GROUPING OF CNAs
grouped_data <- overlap_data_homolog_df %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )

grouped_data <- grouped_data %>%
  group_by(gene_name) %>%
  mutate(count_genes = n()) %>%
  ungroup()

grouped_data <- grouped_data %>%
  group_by(sample_id) %>%
  mutate(count_samples = n()) %>%
  ungroup()

grouped_data <- grouped_data %>%
  group_by(sample_id,gene_name) %>%
  mutate(count_genes_sample = n()) %>%
  ungroup()

grouped_data <- grouped_data %>%
  group_by(gene_name) %>%
  mutate(count_genes_sample_all = sum(count_genes_sample)) %>%
  ungroup()


#| More than 4
grouped_data_f <- grouped_data[which(grouped_data$gene_type == "protein_coding"),]
grouped_data_f <- grouped_data_f[which( (grouped_data_f$count_genes > 3) & (grouped_data_f$count_genes_sample_all == grouped_data_f$count_genes ) ),]

heatmap_plot <- ggplot(grouped_data_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
grouped_data_f$freq_counts <- (grouped_data_f$count_genes/21)*100
grouped_data_f$freq_counts <- as.numeric(grouped_data_f$freq_counts)
grouped_data_f$freq_counts_2 <- round(grouped_data_f$freq_counts)
grouped_data_f$freq_counts_2 <- paste(grouped_data_f$freq_counts_2, "%")
grouped_data_f$freq_counts <- round(grouped_data_f$freq_counts)
grouped_data_f$freq_counts <- as.numeric(grouped_data_f$freq_counts)

#| Unique
unique_gene <- grouped_data_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/HOMOLOG-JOIN-TOP_Candidates-4Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 5.5, width = 8)


#| More than 3
grouped_data_f <- grouped_data[which(grouped_data$gene_type == "protein_coding"),]
grouped_data_f <- grouped_data_f[which( (grouped_data_f$count_genes > 2) & (grouped_data_f$count_genes_sample_all == grouped_data_f$count_genes ) ),]

heatmap_plot <- ggplot(grouped_data_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
grouped_data_f$freq_counts <- (grouped_data_f$count_genes/21)*100
grouped_data_f$freq_counts <- as.numeric(grouped_data_f$freq_counts)
grouped_data_f$freq_counts_2 <- round(grouped_data_f$freq_counts)
grouped_data_f$freq_counts_2 <- paste(grouped_data_f$freq_counts_2, "%")
grouped_data_f$freq_counts <- round(grouped_data_f$freq_counts)
grouped_data_f$freq_counts <- as.numeric(grouped_data_f$freq_counts)

#| Unique
unique_gene <- grouped_data_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/HOMOLOG-JOIN-TOP_Candidates-3Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 7, width = 11)




#########| HEATMAP WITH DATA GROUPING OF CNAs AND CONSIDERING ONLY THE SUB-GROUP

subset <- c("Sample_9", "Sample_8","Sample_4", "Sample_5","Sample_7", "Sample_14", "Sample_23", "Sample_2","Sample_17")
overlap_data_homolog_df_subset <- overlap_data_homolog_df[which(overlap_data_homolog_df$sample_id %in% subset),]
grouped_data_subset <- overlap_data_homolog_df_subset %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )

grouped_data_subset <- grouped_data_subset %>%
  group_by(gene_name) %>%
  mutate(count_genes = n()) %>%
  ungroup()

grouped_data_subset <- grouped_data_subset %>%
  group_by(sample_id) %>%
  mutate(count_samples = n()) %>%
  ungroup()

grouped_data_subset <- grouped_data_subset %>%
  group_by(sample_id,gene_name) %>%
  mutate(count_genes_sample = n()) %>%
  ungroup()

grouped_data_subset <- grouped_data_subset %>%
  group_by(gene_name) %>%
  mutate(count_genes_sample_all = sum(count_genes_sample)) %>%
  ungroup()

#| More than 2
grouped_data_f_subset <- grouped_data_subset[which(grouped_data_subset$gene_type == "protein_coding"),]
grouped_data_f_subset <- grouped_data_f_subset[which( (grouped_data_f_subset$count_genes > 1) & (grouped_data_f_subset$count_genes_sample_all == grouped_data_f_subset$count_genes ) ),]

heatmap_plot <- ggplot(grouped_data_f_subset, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
grouped_data_f_subset$freq_counts <- (grouped_data_f_subset$count_genes/21)*100
grouped_data_f_subset$freq_counts <- as.numeric(grouped_data_f_subset$freq_counts)
grouped_data_f_subset$freq_counts_2 <- round(grouped_data_f_subset$freq_counts)
grouped_data_f_subset$freq_counts_2 <- paste(grouped_data_f_subset$freq_counts_2, "%")
grouped_data_f_subset$freq_counts <- round(grouped_data_f_subset$freq_counts)
grouped_data_f_subset$freq_counts <- as.numeric(grouped_data_f_subset$freq_counts)

#| Unique
unique_gene <- grouped_data_f_subset %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",        # Color for low values (copynumber < 2)
    mid = "white",       # Neutral color (optional, can represent copynumber = 2)
    high = "red",        # Color for high values (copynumber > 2)
    midpoint = 2,        # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/SUBSET-HOMOLOG-JOIN-TOP_Candidates-2Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 5.5, width = 8)



#########| Enrichment 

#| ALL  
overlap_data_homolog_df_f <- overlap_data_homolog_df[which(overlap_data_homolog_df$gene_type == "protein_coding"),]
genes <- unique(overlap_data_homolog_df_f$gene_name[which(overlap_data_homolog_df_f$count_genes >= 1)])
total_gost<- gost(list("Genes with mutations" = genes), 
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
                  custom_bg = custom_bg, 
                  numeric_ns = "", 
                  sources = NULL, 
                  as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste0(dir.results,"Enrichment/ALL_ggplot_results_gostplot_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#writexl::write_xlsx(results,paste("Results/Enrichment/DEL_Enrichment_ggplot_results_gostplot_Protein-Coding_Mut-1Mice_length-", length_cna,"_Exons-", Exons,"_mBIC-",mBIC,".xlsx", sep =""))

dat1_filtered <- results[which(results$source == "TF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/TF_REAC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 9)


dat1_filtered <- results[which(results$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/ALL_REAC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 9)

dat1_filtered <- results[which(results$source == "GO:MF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/ALL_GO-MF_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

dat1_filtered <- results[which(results$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/ALL_GO-CC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 6.5)

dat1_filtered <- results[which(results$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/ALL_GO-BP_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.9)

dat1_filtered <- results[which(results$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/ALL_KEGG_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 7)



#| ALL DELETIONS 
overlap_data_df_f <- overlap_data_df[which(overlap_data_df$gene_type == "protein_coding"),]
overlap_data_df_f <- overlap_data_df_f[which(overlap_data_df_f$type == "del"),]
genes <- unique(overlap_data_df_f$gene_name[which(overlap_data_df_f$count_genes >= 1)])
total_gost<- gost(list("Genes with deletions" = genes), 
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
                  custom_bg = custom_bg, 
                  numeric_ns = "", 
                  sources = NULL, 
                  as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste0(dir.results,"Enrichment/DEL_ggplot_results_gostplot_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#writexl::write_xlsx(results,paste("Results/Enrichment/DEL_Enrichment_ggplot_results_gostplot_Protein-Coding_Mut-1Mice_length-", length_cna,"_Exons-", Exons,"_mBIC-",mBIC,".xlsx", sep =""))

dat1_filtered <- results[which(results$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/DEL_REAC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 9)

dat1_filtered <- results[which(results$source == "GO:MF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/DEL_GO-MF_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

dat1_filtered <- results[which(results$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/DEL_GO-CC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.9)

dat1_filtered <- results[which(results$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/DEL_KEGG_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 7)



#| ALL DUPLICATIONS 
overlap_data_df_f <- overlap_data_df[which(overlap_data_df$gene_type == "protein_coding"),]
overlap_data_df_f <- overlap_data_df_f[which(overlap_data_df_f$type == "dup"),]
genes <- unique(overlap_data_df_f$gene_name[which(overlap_data_df_f$count_genes >= 1)])
total_gost<- gost(list("Genes with Duplications" = genes), 
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
                  custom_bg = custom_bg, 
                  numeric_ns = "", 
                  sources = NULL, 
                  as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste0(dir.results,"Enrichment/DUP_ggplot_results_gostplot_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#writexl::write_xlsx(results,paste("Results/Enrichment/DEL_Enrichment_ggplot_results_gostplot_Protein-Coding_Mut-1Mice_length-", length_cna,"_Exons-", Exons,"_mBIC-",mBIC,".xlsx", sep =""))

dat1_filtered <- results[which(results$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
#ggsave(paste0(dir.results,"Enrichment/DUP_REAC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 9)

dat1_filtered <- results[which(results$source == "GO:MF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/DUP_GO-MF_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 7)

dat1_filtered <- results[which(results$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
#ggsave(paste0(dir.results,"Enrichment/DUP_GO-CC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.9)

dat1_filtered <- results[which(results$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/DUP_GO-BP_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 7.4)




#| ALL HOMOLOG
overlap_data_homolog_df_f <- overlap_data_homolog_df[which(overlap_data_homolog_df$gene_type == "protein_coding"),]
genes <- unique(overlap_data_homolog_df_f$gene_name[which(overlap_data_homolog_df_f$count_genes >= 1)])
total_gost<- gost(list("Genes with mutations" = genes), 
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
                  custom_bg = custom_bg, 
                  numeric_ns = "", 
                  sources = NULL, 
                  as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste0(dir.results,"Enrichment/ALL-HOMOLOG_ggplot_results_gostplot_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#writexl::write_xlsx(results,paste("Results/Enrichment/DEL_Enrichment_ggplot_results_gostplot_Protein-Coding_Mut-1Mice_length-", length_cna,"_Exons-", Exons,"_mBIC-",mBIC,".xlsx", sep =""))

dat1_filtered <- results[which(results$source == "TF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/TF-HOMOLOG_REAC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 9)


dat1_filtered <- results[which(results$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
#ggsave(paste0(dir.results,"Enrichment/ALL-HOMOLOG_REAC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 9)

dat1_filtered <- results[which(results$source == "GO:MF"),]
ggplot(dat1_filtered[1:9,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/ALL-HOMOLOG_GO-MF_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

dat1_filtered <- results[which(results$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/ALL-HOMOLOG_GO-CC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 6.5)

dat1_filtered <- results[which(results$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/ALL-HOMOLOG_GO-BP_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.9)

dat1_filtered <- results[which(results$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
#ggsave(paste0(dir.results,"Enrichment/ALL-HOMOLOG_KEGG_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 7)

dat1_filtered <- results[which(results$source == "HP"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
#ggsave(paste0(dir.results,"Enrichment/ALL-HOMOLOG_KEGG_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 7)


#| ALL DELETIONS 
overlap_data_homolog_df_f <- overlap_data_homolog_df[which(overlap_data_homolog_df$gene_type == "protein_coding"),]
overlap_data_homolog_df_f <- overlap_data_homolog_df_f[which(overlap_data_homolog_df_f$type == "del"),]
genes <- unique(overlap_data_homolog_df_f$gene_name[which(overlap_data_homolog_df_f$count_genes >= 1)])
total_gost<- gost(list("Genes with deletions" = genes), 
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
                  custom_bg = custom_bg, 
                  numeric_ns = "", 
                  sources = NULL, 
                  as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste0(dir.results,"Enrichment/DEL-HOMOLOG_ggplot_results_gostplot_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#writexl::write_xlsx(results,paste("Results/Enrichment/DEL_Enrichment_ggplot_results_gostplot_Protein-Coding_Mut-1Mice_length-", length_cna,"_Exons-", Exons,"_mBIC-",mBIC,".xlsx", sep =""))

dat1_filtered <- results[which(results$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/DEL-HOMOLOG_REAC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 9)

dat1_filtered <- results[which(results$source == "GO:MF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/DEL-HOMOLOG_GO-MF_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

dat1_filtered <- results[which(results$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/DEL-HOMOLOG_GO-CC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.9)

dat1_filtered <- results[which(results$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/DEL-HOMOLOG_KEGG_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 7)



#| ALL DUPLICATIONS 
overlap_data_homolog_df_f <- overlap_data_homolog_df[which(overlap_data_homolog_df$gene_type == "protein_coding"),]
overlap_data_homolog_df_f <- overlap_data_homolog_df_f[which(overlap_data_homolog_df_f$type == "dup"),]
genes <- unique(overlap_data_homolog_df_f$gene_name[which(overlap_data_homolog_df_f$count_genes >= 1)])
total_gost<- gost(list("Genes with Duplications" = genes), 
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
                  custom_bg = custom_bg, 
                  numeric_ns = "", 
                  sources = NULL, 
                  as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste0(dir.results,"Enrichment/DUP-HOMOLOG_ggplot_results_gostplot_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#writexl::write_xlsx(results,paste("Results/Enrichment/DEL_Enrichment_ggplot_results_gostplot_Protein-Coding_Mut-1Mice_length-", length_cna,"_Exons-", Exons,"_mBIC-",mBIC,".xlsx", sep =""))

dat1_filtered <- results[which(results$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
#ggsave(paste0(dir.results,"Enrichment/DUP_REAC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 9)

dat1_filtered <- results[which(results$source == "GO:MF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/DUP-HOMOLOG_GO-MF_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 7)

dat1_filtered <- results[which(results$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
#ggsave(paste0(dir.results,"Enrichment/DUP_GO-CC_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.9)

dat1_filtered <- results[which(results$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
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
ggsave(paste0(dir.results,"Enrichment/DUP_GO-BP_Protein-Coding_Mut-1Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 7.4)




################################################################################


################################################################################
###| CODEX2 metrics evaluation 2
################################################################################

Lkb_cna <- 300
Exons <- 2
mBIC <- 5000
lratio <- 35
raw_cov <- 5
copy_no_del <- 1.7
copy_no_dup <- 2.3

#| Maximun size of the CNA & Ratio between the CNA and the exon
finalcall.CBS_filtered_df <- finalcall.CBS_list_df[which( (finalcall.CBS_list_df$length_kb<=Lkb_cna) & (finalcall.CBS_list_df$length_kb/(finalcall.CBS_list_df$ed_exon-finalcall.CBS_list_df$st_exon+1)<50) ),]

#| Number of exones in the CNA regions |  Minimum lratio
finalcall.CBS_filtered_df <- finalcall.CBS_filtered_df[which(  (finalcall.CBS_filtered_df$lratio>lratio) |  ((finalcall.CBS_filtered_df$ed_exon-finalcall.CBS_filtered_df$st_exon)>=Exons)),]

#| mBIC threshold
finalcall.CBS_filtered_df <- finalcall.CBS_filtered_df[which(finalcall.CBS_filtered_df$mBIC < mBIC),]

#| Copy-number filtering
finalcall.CBS_filtered_df <-finalcall.CBS_filtered_df[which( (finalcall.CBS_filtered_df$copy_no > copy_no_dup) |  (finalcall.CBS_filtered_df$copy_no < copy_no_del) ),]

#| Raw coverage
finalcall.CBS_filtered_df <- finalcall.CBS_filtered_df[which(finalcall.CBS_filtered_df$raw_cov > raw_cov),]


ggplot(finalcall.CBS_filtered_df, aes(x = mBIC, y =norm_cov))+
  geom_point(color="black",size=1)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  ylab("norm_cov")+
  xlab("mBIC")

ggplot(finalcall.CBS_filtered_df, aes(x = mBIC, y =raw_cov))+
  geom_point(color="black",size=1)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  ylab("raw_cov")+
  xlab("mBIC")

length(finalcall.CBS_filtered_df$length_kb[which(finalcall.CBS_filtered_df$length_kb < 1)])
length(finalcall.CBS_filtered_df$length_kb[which(finalcall.CBS_filtered_df$length_kb > 1)])

#| Lkb distribution
ggplot(finalcall.CBS_filtered_df, aes(x = length_kb))+
  geom_histogram(color="black",binwidth=0.1)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("length_kb")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_length_filtered_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

#| lratio distribution
ggplot(finalcall.CBS_filtered_df, aes(x = lratio))+
  geom_histogram(color="black",binwidth=5)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("lratio")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_lratio_filtered_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

#| mBIC distribution
ggplot(finalcall.CBS_filtered_df, aes(x = mBIC))+
  geom_histogram(color="black",binwidth=10)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("mBIC")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_mBIC_filtered_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

#| Raw counts distribution
ggplot(finalcall.CBS_filtered_df, aes(x = raw_cov))+
  geom_histogram(color="black",binwidth=10)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("raw_cov ")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_Raw-Coverage_filtered_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)

#| Normalized raw counts distribution
ggplot(finalcall.CBS_filtered_df, aes(x = norm_cov))+
  geom_histogram(color="black",binwidth=10)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("norm_cov")
ggsave(paste0(dir.results,"CODEX2_Metrics/Histogram_Normalized-Raw-Coverage_filtered_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)




#| Create CNA GRanges object
cna_gr <- GRanges(
  seqnames = finalcall.CBS_filtered_df$chr,        
  ranges = IRanges(start = finalcall.CBS_filtered_df$st_bp, end =finalcall.CBS_filtered_df$ed_bp),
  copynumber = finalcall.CBS_filtered_df$copy_no,
  cnv = finalcall.CBS_filtered_df$cnv,
  length_kb  =finalcall.CBS_filtered_df$length_kb ,
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
  length_kb  = mcols(cna_gr[queryHits(overlaps)])$length_kb ,
  raw_cov = mcols(cna_gr[queryHits(overlaps)])$raw_cov,
  norm_cov = mcols(cna_gr[queryHits(overlaps)])$norm_cov,
  mBIC = mcols(cna_gr[queryHits(overlaps)])$mBIC,
  st_exon = mcols(cna_gr[queryHits(overlaps)])$st_exon,
  ed_exon = mcols(cna_gr[queryHits(overlaps)])$ed_exon
)

#| Filtering places with no gene associated
overlap_data <- overlap_data[which(overlap_data$gene_name != ""),]

#| Filtering genes appearing in the panel of normal
overlap_data <- overlap_data[which( !(overlap_data$gene_name %in%  overlap_data$gene_name[which(overlap_data$sample_id %in% sampname_norm)]) ),]

overlap_data_df <- overlap_data %>%
  group_by(gene_name) %>%
  mutate(count_genes = n()) %>%
  ungroup()

overlap_data_df <- overlap_data_df %>%
  group_by(sample_id,gene_name) %>%
  mutate(count_genes_sample = n()) %>%
  ungroup()

overlap_data_df <- overlap_data_df %>%
  group_by(sample_id) %>%
  mutate(count_samples = n()) %>%
  ungroup()

overlap_data_df <- overlap_data_df %>%
  group_by(gene_name) %>%
  mutate(count_genes_sample_all = sum(count_genes_sample)) %>%
  ungroup()

writexl::write_xlsx(overlap_data_df, paste0(dir.results,"CNVs_results_filtered_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,".xlsx", sep = ""))

gene <- "Pten"
overlap_data_df$type[which(overlap_data_df$gene_name == gene)]
overlap_data_df$sample_id[which(overlap_data_df$gene_name == gene)]

#########| Plotting Histogram
ggplot(overlap_data, aes(x = copynumber, fill = type))+
  geom_histogram(color="black",binwidth=0.1)+
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.text.y = element_text(size = 12, family = "sans")) +
  scale_fill_manual(values = c("blue", "red"))+
  ylab("Counts")+
  xlab("Copy number")
ggsave(paste0(dir.results,"Histograms/Copy_Number_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 5.5)


#########| Plotting barplot
overlap_data_dup_del <- overlap_data %>%
  group_by(sample_id,type) %>%
  summarise(count = n())

ggplot(overlap_data_dup_del, aes(x = reorder(sample_id,-count), y = count, fill = type ))+
  geom_bar(stat = "identity")+
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(color ="black")) +
  xlab("Samples")+
  ylab("Number of genes")+
  scale_fill_manual(values =c( "blue","red"))
ggsave(paste0(dir.results,"Histograms/Dup_Del_Sample_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 6)


############| Chr plot

overlap_data_df$mid_point <- (overlap_data_df$cna_start + overlap_data_df$cna_end)/2

ggplot(overlap_data_df, aes(x = mid_point, y = copynumber , color = copynumber )) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  scale_color_gradient2(low = "blue", 
                        high = "red", 
                        mid = "gray",
                        midpoint = 2) +
  geom_hline(yintercept = 0, alpha =0.5, color ="gray" ) +
  geom_vline(xintercept = seq(1, 100, by = 10), 
             color = "black", linetype = "dotted", alpha = 0.2)+
  geom_hline(yintercept = 2, color ="black") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    strip.text = element_text(size = 8),
    legend.position = "none"
  )  +
  labs(y = "Copy Number", x = "Genomic Position")
ggsave(paste0(dir.results,"Heatmaps/CNAs_Chromosome_Position-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 12)


#| Plotting CNVs by chromosomes
unique_cna <- overlap_data_df %>% distinct(cna_index, .keep_all = TRUE)

new_df <- unique_cna[,c("chr","cna_start", "cna_end", "copynumber", "sample_id")]
names(new_df) <- c("chromosome", "start", "end", "segmean", "sample")
new_df$chromosome <- gsub("chr", "", new_df$chromosome)
new_df$segmean <- as.integer(new_df$segmean)

cnFreq(new_df, genome="mm10",
       CN_Loss_colour="blue",
       CN_Gain_colour ="red",
       CN_low_cutoff = 1.7,
       CN_high_cutoff = 2.3) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        text=element_text(size=9,  family="sans"),
        axis.title.y = element_text(face = "bold", size =12),
        axis.title.x = element_text(size =12))
ggsave(paste0(dir.results,"Heatmaps/cnFreq_Proportion_CopyNumber_Gain_Loss-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 4.5, width = 15)

cnSpec(
  x = new_df, 
  genome = "mm10", 
  plot_title = "Global Copy Number Alterations (CNA)",
  CN_Loss_colour = "blue", 
  CN_Gain_colour = "red"
)+theme_minimal() +  # Apply a minimal theme
  theme(
    legend.title = element_blank(),  # Remove legend title
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    strip.text.y=element_text(angle=0)
  )
ggsave(paste0(dir.results,"Heatmaps/cnSpec_Gain_Loss-Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 5.5, width = 9)


#########|  HEATMAP top candidates More than 6
overlap_data_df_f <- overlap_data_df[which(overlap_data_df$gene_type == "protein_coding"),]
overlap_data_df_f <- overlap_data_df_f[which( (overlap_data_df_f$count_genes > 5) & (overlap_data_df_f$count_genes_sample_all == overlap_data_df_f$count_genes ) ),]

overlap_data_df_f <- overlap_data_df_f %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon,count_genes,count_genes_sample_all) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )

heatmap_plot <- ggplot(overlap_data_df_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
overlap_data_df_f$freq_counts <- (overlap_data_df_f$count_genes/21)*100
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- paste(overlap_data_df_f$freq_counts_2, "%")
overlap_data_df_f$freq_counts <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)

#| Unique
unique_gene <- overlap_data_df_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/TOP_Candidates-6Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 6.5, width = 10)


#########|  HEATMAP top candidates More than 5
overlap_data_df_f <- overlap_data_df[which(overlap_data_df$gene_type == "protein_coding"),]
overlap_data_df_f <- overlap_data_df_f[which( (overlap_data_df_f$count_genes > 4) & (overlap_data_df_f$count_genes_sample_all == overlap_data_df_f$count_genes ) ),]

overlap_data_df_f <- overlap_data_df_f %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon,count_genes,count_genes_sample_all) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )

heatmap_plot <- ggplot(overlap_data_df_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
overlap_data_df_f$freq_counts <- (overlap_data_df_f$count_genes/21)*100
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- paste(overlap_data_df_f$freq_counts_2, "%")
overlap_data_df_f$freq_counts <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)

#| Unique
unique_gene <- overlap_data_df_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/TOP_Candidates-5Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 8, width = 10)


#########|  HEATMAP top candidates More than 4
overlap_data_df_f <- overlap_data_df[which(overlap_data_df$gene_type == "protein_coding"),]
overlap_data_df_f <- overlap_data_df_f[which( (overlap_data_df_f$count_genes > 3) & (overlap_data_df_f$count_genes_sample_all == overlap_data_df_f$count_genes ) ),]
overlap_data_df_f <- overlap_data_df_f %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon,count_genes,count_genes_sample_all) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )

heatmap_plot <- ggplot(overlap_data_df_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
overlap_data_df_f$freq_counts <- (overlap_data_df_f$count_genes/21)*100
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- paste(overlap_data_df_f$freq_counts_2, "%")
overlap_data_df_f$freq_counts <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)

#| Unique
unique_gene <- overlap_data_df_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/TOP_Candidates-4Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 14, width = 12)


#########|  HEATMAP top candidates More than 3
overlap_data_df_f <- overlap_data_df[which(overlap_data_df$gene_type == "protein_coding"),]
overlap_data_df_f <- overlap_data_df_f[which( (overlap_data_df_f$count_genes > 2) & (overlap_data_df_f$count_genes_sample_all == overlap_data_df_f$count_genes ) ),]
overlap_data_df_f <- overlap_data_df_f %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon,count_genes,count_genes_sample_all) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )

heatmap_plot <- ggplot(overlap_data_df_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
overlap_data_df_f$freq_counts <- (overlap_data_df_f$count_genes/21)*100
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts_2 <- paste(overlap_data_df_f$freq_counts_2, "%")
overlap_data_df_f$freq_counts <- round(overlap_data_df_f$freq_counts)
overlap_data_df_f$freq_counts <- as.numeric(overlap_data_df_f$freq_counts)

#| Unique
unique_gene <- overlap_data_df_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/TOP_Candidates-3Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 25, width = 13)


#########| HEATMAP WITH ONLY HUMAN HOMOLOG
homolog_info <- getBM(
  attributes = c("external_gene_name",  # Mouse Ensembl ID (for merging)
                 "hsapiens_homolog_associated_gene_name", 
                 "hsapiens_homolog_ensembl_gene", 
                 "hsapiens_homolog_orthology_type"),
  mart = ensembl
)
names(homolog_info) <- c(c("gene_name",  # Mouse Ensembl ID (for merging)
                           "hsapiens_homolog_associated_gene_name", 
                           "hsapiens_homolog_ensembl_gene", 
                           "hsapiens_homolog_orthology_type"))
annotations_homolog <- merge(annotations, homolog_info, by = "gene_name")
annotations_homolog <- annotations_homolog[which(annotations_homolog$hsapiens_homolog_associated_gene_name !=""),]

#| Selecting only homologs
overlap_data_homolog_df <- overlap_data_df[which(overlap_data_df$gene_name %in% unique(annotations_homolog$gene_name)),]


#########|  HEATMAP top candidates More than 4
overlap_data_homolog_df_f <- overlap_data_homolog_df[which(overlap_data_homolog_df$gene_type == "protein_coding"),]
overlap_data_homolog_df_f <- overlap_data_homolog_df_f[which( (overlap_data_homolog_df_f$count_genes > 3) & (overlap_data_homolog_df_f$count_genes_sample_all == overlap_data_homolog_df_f$count_genes ) ),]
overlap_data_homolog_df_f <- overlap_data_homolog_df_f %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon,count_genes,count_genes_sample_all) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )

heatmap_plot <- ggplot(overlap_data_homolog_df_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
overlap_data_homolog_df_f$freq_counts <- (overlap_data_homolog_df_f$count_genes/21)*100
overlap_data_homolog_df_f$freq_counts <- as.numeric(overlap_data_homolog_df_f$freq_counts)
overlap_data_homolog_df_f$freq_counts_2 <- round(overlap_data_homolog_df_f$freq_counts)
overlap_data_homolog_df_f$freq_counts_2 <- paste(overlap_data_homolog_df_f$freq_counts_2, "%")
overlap_data_homolog_df_f$freq_counts <- round(overlap_data_homolog_df_f$freq_counts)
overlap_data_homolog_df_f$freq_counts <- as.numeric(overlap_data_homolog_df_f$freq_counts)

#| Unique
unique_gene <- overlap_data_homolog_df_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/HOMOLOG-TOP_Candidates-4Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 13, width = 12)


#########|  HEATMAP top candidates More than 3
overlap_data_homolog_df_f <- overlap_data_homolog_df[which(overlap_data_homolog_df$gene_type == "protein_coding"),]
overlap_data_homolog_df_f <- overlap_data_homolog_df_f[which( (overlap_data_homolog_df_f$count_genes > 2) & (overlap_data_homolog_df_f$count_genes_sample_all == overlap_data_homolog_df_f$count_genes ) ),]
overlap_data_homolog_df_f <- overlap_data_homolog_df_f %>%
  group_by(chr, cna_index, cna_width, cna_start, cna_end,gene_type, copynumber, type, lratio, sample_id, length_kb, raw_cov, norm_cov, mBIC, st_exon, ed_exon,count_genes,count_genes_sample_all) %>%
  summarise(
    gene_name = paste(unique(gene_name), collapse = ", "),
    .groups = "drop"
  )
heatmap_plot <- ggplot(overlap_data_homolog_df_f, aes(x = reorder(sample_id,-count_genes), y = reorder(gene_name, count_genes), fill= copynumber)) + 
  geom_tile(color ="grey40") +
  theme_minimal() +
  theme(text=element_text(size=14,  family="sans"), 
        axis.text.x = element_text(size=10,color ="black",angle = 60, hjust=1),
        axis.text.y = element_text(size=10,color ="black"), 
        legend.position = "none") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  )+
  xlab("Samples")+
  ylab("Genes") +
  labs(fill = "Call") 

#| Computing the frequency
overlap_data_homolog_df_f$freq_counts <- (overlap_data_homolog_df_f$count_genes/21)*100
overlap_data_homolog_df_f$freq_counts <- as.numeric(overlap_data_homolog_df_f$freq_counts)
overlap_data_homolog_df_f$freq_counts_2 <- round(overlap_data_homolog_df_f$freq_counts)
overlap_data_homolog_df_f$freq_counts_2 <- paste(overlap_data_homolog_df_f$freq_counts_2, "%")
overlap_data_homolog_df_f$freq_counts <- round(overlap_data_homolog_df_f$freq_counts)
overlap_data_homolog_df_f$freq_counts <- as.numeric(overlap_data_homolog_df_f$freq_counts)

#| Unique
unique_gene <- overlap_data_homolog_df_f %>% distinct(sample_id,gene_name, .keep_all = TRUE)
unique_gene$count_genes_2 <- 1
unique_gene <- as.data.frame(unique_gene)

#| Barplot
barplot_plot <- ggplot(unique_gene, aes(x = count_genes_2 , y = reorder(gene_name, freq_counts), fill = copynumber)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Percentage") +
  scale_fill_gradient2(
    low = "blue",       # Color for low values (copynumber < 2)
    mid = "white",      # Neutral color (optional, can represent copynumber = 2)
    high = "red",       # Color for high values (copynumber > 2)
    midpoint = 2,       # Define the midpoint for the gradient
    name = "Copy Number" # Legend title
  ) +
  ylab("Number of Genes") +
  labs(fill = "CopyNumber") +
  geom_text(aes(x = count_genes, label = freq_counts_2),hjust = -0.5, size = 2, color = "black") +
  xlim(0,30)

heatmap_plot + barplot_plot + plot_layout(ncol = 2,widths = c(4, 1))
ggsave(paste0(dir.results,"Heatmaps/HOMOLOG-TOP_Candidates-3Mice_Lkb-", Lkb_cna,"_Exons-", Exons,"_mBIC-",mBIC,"_copy_no-Del-",copy_no_del,"-Amp-", copy_no_dup,".pdf"), height = 23, width = 12)

