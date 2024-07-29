#!/bin/bash
#$ -N scheduler
#$ -cwd
#$ -j y


micromamba activate varrna
snakemake -p \
    --cluster "qsub -pe smp {threads} -cwd -o qsub_logfiles -e qsub_logfiles" \
    --jobname "sm_{rule}_{wildcards.sample}_{jobid}" \
    --jobs 64 \
    --rerun-incomplete \
    --latency-wait 10 
