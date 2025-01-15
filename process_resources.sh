#!/bin/bash

set -euo pipefail
trap 'echo "Error: Script failed at line $LINENO."; exit 1' ERR

resources="${1:-resources}"

# Create dictionary and index for reference genome
gatk CreateSequenceDictionary -R $resources/GRCh38.p13.genome.fa
samtools faidx $resources/GRCh38.p13.genome.fa

# Salmon index
dependencies/salmon-1.9.0_linux_x86_64/bin/salmon index -t $resources/gencode.v47.transcripts.fa.gz -i dependencies/gencode.v47.transcripts_index

# Gene names bed file
# export some of the gtf information to BED file
awk 'BEGIN{OFS="\t"} $3 == "gene" {
    gene_name = "UNDEFINED";
    for (i = 9; i <= NF; i++) {
        if ($i == "gene_name") {
            gname_field = $(i+1);
            gsub("\"", "", gname_field);
            gsub(";", "", gname_field);
            gene_name = gname_field;
            print $1, $4, $5, gene_name;
            break;
        }
    }
}' $resources/gencode.v43.primary_assembly.annotation.gtf > $resources/gencode.v43.primary_assembly.gene_name.bed

## merge gene regions of multiple entries
sort -k1,1 -k4,4 $resources/gencode.v43.primary_assembly.gene_name.bed | \
bedtools groupby -i stdin -g 1,4 -c 2,3 -o min,max | \
    awk -v OFS="\t" '{print $1,$3,$4,$2}' | \
    sort -k1,1 -k2,2n > \
    $resources/tmp && mv $resources/tmp $resources/gencode.v43.primary_assembly.gene_name.bed

bgzip $resources/gencode.v43.primary_assembly.gene_name.bed
tabix -p bed $resources/gencode.v43.primary_assembly.gene_name.bed.gz

# dbSNP known sites for GATK BaseRecalibrator
zcat $resources/00-common_all.vcf.gz | sed '/^[0-9XY]/s/^/chr/' | bgzip > $resources/dbsnp151_common.hg38.vcf.gz
tabix -p vcf $resources/dbsnp151_common.hg38.vcf.gz
rm $resources/00-common_all.vcf.gz

# Repeat Masker
zcat $resources/hg38.fa.out.gz | awk '{OFS="\t"} { if (NR>3) {print $5, $6, $7, $11, "0"} }' > $resources/repmask_hg38.bed
rm $resources/hg38.fa.out.gz

# REDI Portal
cat $resources/vcf_header.txt > $resources/RNAedit.vcf
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> $resources/RNAedit.vcf
zcat $resources/TABLE1_hg38_v3.txt.gz | awk '{OFS="\t"} { if (NR>1) {print $2, $3, ".", "A", "G", ".", ".", "."} }' >> $resources/RNAedit.vcf

bgzip $resources/RNAedit.vcf
bcftools sort -Oz -o $resources/RNAedit.sorted.vcf.gz $resources/RNAedit.vcf.gz
tabix -p vcf $resources/RNAedit.sorted.vcf.gz
rm $resources/RNAedit.vcf.gz $resources/TABLE1_hg38_v3.txt.gz