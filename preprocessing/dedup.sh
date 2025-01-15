#!/bin/bash

# call the script together with the bamfile

bam=$1

fumi_tools dedup --paired --threads 8 --memory 100G -i $bam -o ${bam%.bam}.dedup.bam 

samtools sort --threads 8 -O BAM -o ${bam%.bam}.dedup.srt.bam ${bam%.bam}.dedup.bam
samtools index -@ 8 ${bam%.bam}.dedup.srt.bam

rm ${bam%.bam}.dedup.bam
