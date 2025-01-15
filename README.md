# sRNAseq_varmiR
In this repository, I include the scripts used to preprocess small RNA sequencing (sRNAseq) data. In addition, I include scripts developed to detect miRNA with variable expression patterns in a group versus another

## Map

```
sbatch -t 05:00:00 --mem=50G --wrap="./fastMGcount/submit_mgcount.sh T2DvsND_MGcount T2DvsND list_bamfiles.txt ./mgcount_human/Homo_sapiens.GRCh38.custom.gtf"
```
