#!/bin/bash

fasta="../resources/GRCh38.p13.genome.fa"
gtf="../resources/gencode.v43.primary_assembly.annotation.gtf"


# sjdbOverhang should be max(read length)-1
# https://www.biostars.org/p/93883/
dependencies/STAR --runMode genomeGenerate \
    --runThreadN 20 \
    --genomeDir ../resources/star_index \
    --genomeFastaFiles $fasta \
    --sjdbGTFfile $gtf \
    --sjdbOverhang 160



