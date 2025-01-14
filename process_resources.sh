#!/bin/bash

resources="${1:-resources}"

# Gene names bed file
## export some of the gtf information to BED file
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
zcat $resources/rmsk.txt.gz | awk '{OFS="\t"} { {print $6, $7, $8, $12, "0", $10} }' > $resources/repmask_hg38.bed

# REDI Portal
cat $resources/vcf_header.txt > $resources/RNAedit.vcf
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> $resources/RNAedit.vcf

zcat resources/TABLE1_hg38.txt.gz | awk '{OFS="\t"} { if (NR>1) {print $1, $2, ".", "A", "G", ".", ".", "."} }' >> resources/RNAedit.vcf

bgzip $resources/RNAedit.vcf
bcftools sort -Oz -o $resources/RNAedit.sorted.vcf.gz $resources/RNAedit.vcf.gz
tabix -p vcf $resources/RNAedit.sorted.vcf.gz
rm $resources/RNAedit.vcf.gz