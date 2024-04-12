#!/bin/bash

SAMPLE=$1
VCFDIR=$2

vcf_files=()
for i in {1..22} {X,Y}
do
    vcf_files+=("$VCFDIR/$SAMPLE.chr$i.vcf.gz")
done

bcftools concat --threads 20 --naive -Oz -o $VCFDIR/$SAMPLE.vcf.gz "${vcf_files[@]}"
tabix -p vcf $VCFDIR/$SAMPLE.vcf.gz