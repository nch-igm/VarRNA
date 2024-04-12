#!/bin/bash

# Intersect vcf with coverage bed file

# Filter VCF by 
## - no repetitive elements
## - exome regions
## - biallelic PASS or '.' variants only

INPUT=$1
OUTPUT=$2
BEDREP=$3
EXOME=$4

bedtools intersect -a $INPUT -b $BEDREP -v -header | \
    bedtools intersect -a - -b $EXOME -header | \
    bcftools view --threads 8 -f.,PASS -i 'MIN(FMT/DP)>10 & QUAL>100 & QD>2' -m2 -M2 -Oz -o $OUTPUT -

tabix -p vcf $OUTPUT

