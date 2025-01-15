#!/bin/bash

outdir=$1
outname=$2
inbamfile=$3

path2featureCounts=~/bin/subread-2.0.1-source/bin/
path2samtools=~/bin/samtools-1.11/

while read -r line
do
    inbam=$(echo $line | awk '{print $1}')
    sample=${outdir}/$(echo $line | awk '{print $2}')
    outtmp=$(echo $inbam | awk '{for (i=length($1); i>=1; i--) {if (substr($1,i,1)=="/") print substr($1,i+1,length($1))}}' | head -1)
    outtmp=$(echo $outtmp | awk -v outtmp=${outtmp} -v inbam=${inbam} '{if (length(outtmp)==0) {print inbam} else {print outtmp}}')
    echo $sample $outtmp
    ############
    # Counting #
    ############

    # Round 1: counting small non coding hits
    ${path2featureCounts}/featureCounts -o ${sample}_small.csv -a ${outdir}/${outname}_small_annot.gtf -F GTF -t transcript -g transcript_name -O --fracOverlap 1 --largestOverlap -M --fraction -s 1 -R BAM $inbam
    mv ${outdir}/${outtmp}.featureCounts.bam ${sample}.small.bam
    ${path2samtools}/samtools view -h ${sample}.small.bam | grep -v 'XS:Z:Assigned' | samtools view -O BAM -o ${sample}.small_notAssigned.bam

    # Round 2: counting long non coding exonic hits
    ${path2featureCounts}/featureCounts -o ${sample}_long_exons.csv -a ${outdir}/${outname}_long_annot.gtf -F GTF -t exon -g gene_name -O --fracOverlap 1 -M --fraction -s 1 -R BAM ${sample}.small_notAssigned.bam
    mv ${sample}.small_notAssigned.bam.featureCounts.bam ${sample}.longexon.bam
    ${path2samtools}/samtools view -h ${sample}.longexon.bam | grep -v 'XS:Z:Assigned' | samtools view -O BAM -o ${sample}.longexon_notAssigned.bam

    # Round 3: counting long non coding intronic hits
    ${path2featureCounts}/featureCounts -o ${sample}_long_introns.csv -a ${outdir}/${outname}_long_annot.gtf -F GTF -t gene -g gene_name -O --fracOverlap 1 -M --fraction -s 1 ${sample}.longexon_notAssigned.bam

    # for later:
    ${path2samtools}/samtools view ${sample}.small.bam | awk '/XS:Z:Assigned/ {print $1"\t"substr($NF,6,length($NF))}' > ${sample}.fc_small.tsv

    # clean up:
    rm ${sample}.small.bam ${sample}.small_notAssigned.bam ${sample}.longexon.bam ${sample}.longexon_notAssigned.bam

done < $inbamfile
