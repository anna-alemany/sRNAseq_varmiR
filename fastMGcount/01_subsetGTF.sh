#!/bin/bash

outdir=$1
outname=$2 # root for output table
gtf=$3 # ./fastMGcount/Homo_sapiens.GRCh38.custom.gtf

mkdir -p $outdir

small='transcript_biotype "snRNA"|transcript_biotype "snoRNA"|transcript_biotype "scaRNA"|transcript_biotype "pre_miRNA"|transcript_biotype "miRNA"|transcript_biotype "Y_RNA"|transcript_biotype "siRNA"|transcript_biotype "vault_RNA"|transcript_biotype "ribozyme"|transcript_biotype "misc_RNA"|transcript_biotype "piRNA"|transcript_biotype "scRNA"|transcript_biotype "sRNA"|transcript_biotype "rRNA"|transcript_biotype "rRNA_pseudogene"|transcript_biotype "Mt_rRNA"|transcript_biotype "tRNA"|transcript_biotype "tRF"|transcript_biotype "tRF5"|transcript_biotype "tRF3"|transcript_biotype "Mt_tRNA"'

awk -v pat="$small" '$0 ~ pat' $gtf > ${outdir}/${outname}_small_annot.gtf

long='gene_biotype "ERCC"|gene_biotype "IG_C_gene"|gene_biotype "IG_D_gene"|gene_biotype "IG_J_gene"|gene_biotype "IG_V_gene"|gene_biotype "IG_LV_gene"|gene_biotype "IG_pseudogene"|gene_biotype "IG_C_pseudogene"|gene_biotype "IG_D_pseudogene"|gene_biotype "IG_J_pseudogene"|gene_biotype "IG_V_pseudogene"|gene_biotype "TR_C_gene"|gene_biotype "TR_D_gene"|gene_biotype "TR_J_gene"|gene_biotype "TR_V_gene"|gene_biotype "TR_J_pseudogene"|gene_biotype "TR_V_pseudogene"|gene_biotype "protein_coding"|gene_biotype "lncRNA"|gene_biotype "lincRNA"|gene_biotype "macro_lncRNA"|gene_biotype "ncRNA"|gene_biotype "antisense"|gene_biotype "sense_intronic"|gene_biotype "sense_overlapping"|gene_biotype "processed_transcript"|gene_biotype "3prime_overlapping_ncRNA"|gene_biotype "bidirectional_promoter_lncRNA"|gene_biotype "nontranslating_CDS"|gene_biotype "antisense_RNA"|gene_biotype "pseudogene"|gene_biotype "polymorphic_pseudogene"|gene_biotype "unitary_pseudogene"|gene_biotype "unprocessed_pseudogene"|gene_biotype "processed_pseudogene"|gene_biotype "transcribed_unprocessed_pseudogene"|gene_biotype "transcribed_processed_pseudogene"|gene_biotype "transcribed_unitary_pseudogene"|gene_biotype "translated_unprocessed_pseudogene"|gene_biotype "translated_processed_pseudogene"|gene_biotype "nontranslating_CDS"'

awk -v pat="$long" '$0 ~ pat' $gtf > ${outdir}/${outname}_long_annot.gtf
