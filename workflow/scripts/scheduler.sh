#!/bin/bash
#$ -N scheduler
#$ -cwd
#$ -j y


conda activate tetra
snakemake -p \
    --cluster "qsub -pe smp {threads} -cwd -o qsub_logfiles -e qsub_logfiles" \
    --jobname "sm_{rule}_{wildcards.sample}_{jobid}" \
    -j 48 \
    --rerun-incomplete \
    --latency-wait 60 
