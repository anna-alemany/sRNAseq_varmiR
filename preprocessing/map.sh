#!/bin/bash

# provide input fastq file to map.
# It is assumed that STAR is in your bin folder

inputfq=$1

genome=./ensembl/human/99/star_v277a_primaryAssemblyERCC_index_150 ### You might have to modify this according to your needs!
outprefix=${inputfq%.fastq.gz}

STAR --runThreadN 8 --genomeDir ${genome} --readFilesIn ${inputfq} --readFilesCommand zcat --outFilterMultimapNmax 20 --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMattributes All  --outFileNamePrefix ${outprefix}

