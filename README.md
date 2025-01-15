# sRNAseq_varmiR
In this repository, I include the scripts used to preprocess small RNA sequencing (sRNAseq) data. In addition, I include scripts developed to detect miRNA with variable expression patterns in a group versus another

## Fastq file preprocessing, mapping and deduplication

We assume that single end sequencing was performed. Our reads start with a 12nt long UMI, according to the standard practices in Diagenode (https://www.diagenode.com/en/p/D-Plex-Small-RNA-seq-Library-Prep-x24).

1) Trimming
2) Removing UMI sequence to the read, and appending it to the read name with a "_" (for posterior use of umitools)
3) Ribosomal depletion (`./preprocessing/ribodep.sh $fastq`)

4) Mapping (`./preprocessing/map.sh $fastq`)
   This has been done using STAR. We recommend the genome and gtf used in the MGcount original paper (https://filedn.com/lTnUWxFTA93JTyX3Hvbdn2h/mgcount/UserGuide.html), which can be obtained by:
   ```
   wget https://filedn.com/lTnUWxFTA93JTyX3Hvbdn2h/mgcount/integrated_annotations_gtf.zip
   unzip integrated_annotations_gtf.zip -d annotations_gtf
   ```

6) Deduplication (`./preprocessing/dedup.sh $bamfile`)
   This is done with fumi_tools (https://ccb-gitlab.cs.uni-saarland.de/tobias/fumi_tools)

## Count table with fastMGcount

For this project, we wrote a faster version of MGcount, that we call fastMGcount. It can be found in this repository. We still use the GTF that was published with the original software (see above). We don't provide the GTF in this repository for size reasons. 
To call the fastMGcount pipeline, one can use the provided `submit_mgcount.sh` wrapper as follows: 

```
sbatch -t 05:00:00 --mem=50G --wrap="./fastMGcount/submit_mgcount.sh ${output_directory} ${output_name} ${list_bamfiles} ${gtf_file}"
```

This are the variables:
* ${output_directory}: string with name for output directory
* ${output_name}: string providing prefix for output files
* ${list_bamfiles}: file containing all BAM files and short names to be provided to the fastMGcount software. It is a tsv file, that reads as:
   ```
   R155ND_trimmed_R1.nonRiboAligned.sortedByCoord.out.dedup.srt.bam	R155ND
   R167ND_trimmed_R1.nonRiboAligned.sortedByCoord.out.dedup.srt.bam	R167ND
   R171ND_trimmed_R1.nonRiboAligned.sortedByCoord.out.dedup.srt.bam	R171ND
   R185T2D_trimmed_R1.nonRiboAligned.sortedByCoord.out.dedup.srt.bam	R185T2D
   R195T2D_trimmed_R1.nonRiboAligned.sortedByCoord.out.dedup.srt.bam	R195T2D
   R201T2D_trimmed_R1.nonRiboAligned.sortedByCoord.out.dedup.srt.bam	R201T2D
   ```
* {gtf_file}: gtf file (see above)

## Identification of variable miRNA in T2D
Using DESeq2 (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8), we obtained median-of-ratios normalized count tables. This can be achieved using the script: 
```
./var_miRNA/deseq2.R
```
, that uses the count tables provided by fastMGcount. 

Next, we performed Principal Component analysis and calculated variance and corresponding error bars in miRNA expression per group to identify variable miRNA in each group. This can be achieved by running the script: 
```
./var_miRNA/PCA_diff_expr.py
```
