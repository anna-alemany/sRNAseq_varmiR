#!/bin/bash

outdir=$1
outname=$2
inbamfile=$3
gtf=$4

./fastMGcount/01_subsetGTF.sh ${outdir} ${outname} ${gtf}

./fastMGcount/02_readAssignation.sh ${outdir} ${outname} ${inbamfile}

for sample in $(awk '{print $NF}' $inbamfile)
do
    ./fastMGcount/03_mg_small.py ${outdir}/${sample} ${outdir}/${outname}
done

./fastMGcount/04_communities_small.py ${outdir} ${outname} ${inbamfile}
