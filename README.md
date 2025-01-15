# sRNAseq_varmiR
In this repository, I include the scripts used to preprocess small RNA sequencing (sRNAseq) data. In addition, I include scripts developed to detect miRNA with variable expression patterns in a group versus another

## Fastq file preprocessing, mapping and deduplication

We assume that single end sequencing was performed. Our reads start with a 12nt long UMI, according to the standard practices in Diagenode (https://www.diagenode.com/en/p/D-Plex-Small-RNA-seq-Library-Prep-x24).

1) Trimming
2) Removing UMI sequence to the read, and appending it to the read name with a "_" (for posterior use of umitools)
3) Ribosomal depletion (`./preprocessing/ribodep.sh $fastq`)

4) Mapping

5) Deduplication (`./preprocessing/dedup.sh $bamfile`)
   This is done with fumi_tools (https://ccb-gitlab.cs.uni-saarland.de/tobias/fumi_tools)

## Count table with fastMGcount

```
sbatch -t 05:00:00 --mem=50G --wrap="./fastMGcount/submit_mgcount.sh T2DvsND_MGcount T2DvsND list_bamfiles.txt ./mgcount_human/Homo_sapiens.GRCh38.custom.gtf"
```
