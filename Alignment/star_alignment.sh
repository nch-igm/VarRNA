#!/bin/bash

fasta="../resources/GRCh38.p13.genome.fa"
gtf="../resources/gencode.v43.primary_assembly.annotation.gtf"

input_dir=../../RNA_VC/public_data_for_varrna_example
r1=$input_dir/${1}_1.fastq.gz
r2=$input_dir/${1}_2.fastq.gz

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

