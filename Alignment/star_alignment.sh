#!/bin/bash

fasta="../resources/GRCh38.p13.genome.fa"
gtf="../resources/gencode.v43.primary_assembly.annotation.gtf"

r1=/igm/projects/Audrey_RNAseq_Analysis/data/nch-igm-research-cancer/$1/RNASeq/FASTQ/$1.R1.fastq.gz
r2=/igm/projects/Audrey_RNAseq_Analysis/data/nch-igm-research-cancer/$1/RNASeq/FASTQ/$1.R2.fastq.gz

outdir=../results/$1/BAMs


# First pass
dependencies/STAR \
    --runThreadN 20 \
    --genomeDir ../resources/star_index \
    --readFilesIn $r1 $r2 \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $outdir/$1.temp. \
    --outSAMunmapped Within

# Second pass
dependencies/STAR \
    --runThreadN 20 \
    --genomeDir ../resources/star_index \
    --readFilesIn $r1 $r2 \
    --readFilesCommand zcat \
    --sjdbFileChrStartEnd $outdir/$1.temp.SJ.out.tab \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $outdir/$1. \
    --outSAMunmapped Within

