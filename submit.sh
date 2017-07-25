#! /bin/bash

mkdir -p $(pwd)/_log

snakemake -pn \
    --cluster 'qsub -v PATH="/home/miniconda3/bin:$PATH" -j y -o $(pwd)/_log/ -b n -l {params.req} -S /bin/bash' \
    --jobs 99 \
    --latency-wait 60 \
    --restart-times 0 \
    --jn "sk.{jobid}.{rulename}.sh" \
    --use-conda
