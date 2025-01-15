#!/bin/bash
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# bwa and samtools are assumed to be in the user's path.                 #
# USAGE: Provide fastq file name as an input and let the script work!    #
# Submission with 8 threads by default is assumed                        #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

fastq=$1 # R155ND_trimmed_R1.fastq.gz
out=${fastq%.f*q.gz}

riboref=./preprocessing/rRNA.homo_sapiens.fa

bwa aln $riboref $fastq > ${out}.aln.sai
bwa samse $riboref ${out}.aln.sai $fastq | samtools view -Sb > ${out}.aln-ribo.bam &

bwa mem -t 8 -h 15 $riboref $fastq | samtools view -Sb > ${out}.mem-ribo.bam &

wait

samtools merge -n -r -h ${out}.aln-ribo.bam --threads 8 ${out}.all-ribo.bam ${out}.aln-ribo.bam ${out}.mem-ribo.bam
samtools sort -n --threads 8 ${out}.all-ribo.bam -O BAM -o ${out}.nsorted.all-ribo.bam

samtools view ${out}.nsorted.all-ribo.bam | awk '
        BEGIN {OFS="\n"; name="init"; count=1; flagsum=0; read = "init"; seq = "NNN"; phred = "###"} {
            if (name != $1) {if (flagsum == 4*count) {print "@"name, seq, "+", phred}; name=$1; count = 0; flagsum = 0; seq = $10; phred = $11};
            count += 1; if ($2 != 0) {flagsum += 4}
        } END {
            if (flagsum == 4*count) {print "@"name, seq, "+", phred}
        }' > ${out}.nonRibo.fastq

gzip ${out}.nonRibo.fastq

rm ${out}.aln-ribo.bam
rm  ${out}.mem-ribo.bam
rm ${out}.all-ribo.bam
rm ${out}.aln.sai
rm ${out}.nsorted.all-ribo.bam
