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

```
sbatch -t 05:00:00 --mem=50G --wrap="./fastMGcount/submit_mgcount.sh T2DvsND_MGcount T2DvsND list_bamfiles.txt ./fastMGcount/Homo_sapiens.GRCh38.custom.gtf"
```
