"deseq2.R" 44L, 1666C                                                                                                                                                1,1           All
#!/usr/bin/Rscript
# to run this: Rscript deseq2.R (keep an eye on paths)

library(dplyr)

# Input data
srna_metadata_df <- read.csv('./output_fastMGcount/sample_metadata.csv', sep = ',') # output file provided by fastMGcount
srna_metadata_df$type <- relevel(factor(srna_metadata_df$type), ref = 'ND')
srna_matrix <- read.csv('./output_fastMGcount/merged_miRNA_aggregated_counts.tsv', sep = '\t', row.names = 1) # output file provided by fastqMGcount

rownames(srna_metadata_df) <- srna_metadata_df$sampleID
colnames(srna_matrix) <- rownames(srna_metadata_df)

# feature selection
filter_lowfeatures <- function(M, filt_ns, filt_th, cpm = TRUE){
  nM <- M; if(cpm) nM <- apply(M, 2, function(x) x/sum(x)*10^6)
  inc <- apply(nM, 1, function(x) sum(as.numeric(x)> filt_th) >= filt_ns)
  return(inc)}

thS <- 3
th_cpm_small <- 3
incl_f <- filter_lowfeatures(round(srna_matrix), filt_ns = thS, filt_th = th_cpm_small, cpm = TRUE)
srna_matrix <- srna_matrix[incl_f,]

# deseq2
srna_dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(srna_matrix),
  colData = srna_metadata_df,
  design = ~ type)

srna_dds <- DESeq2::DESeq(srna_dds, parallel = TRUE)

srna_vst1 <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(srna_dds, blind = FALSE))
srna_vst2 <- t(apply(srna_vst1, 1, function(x) tapply(x, srna_metadata_df$type, mean)))
srna_norm <- DESeq2::counts(srna_dds, normalized=TRUE)

res <- DESeq2::results(srna_dds)

save.image(file = 'mirna_dds.Rdata')

write.table(srna_vst1, file = 'deseq2_vst_data.tsv', sep = '\t')
write.table(srna_norm, file = 'deseq2_norm_data.tsv', sep = '\t')
write.table(res, file = 'deseq_dextable.tsv', sep = '\t')
