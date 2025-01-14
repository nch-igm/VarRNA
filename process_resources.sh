#!/bin/bash

resources="${1:-resources}"

####### Create salmon index for transcript quantification ########
dependencies/salmon-1.9.0_linux_x86_64/bin/salmon index -t $resources/gencode.v47.transcripts.fa.gz -i dependencies/gencode.v47.transcripts_index

####### Export some of the gtf information to BED file and merge gene regions of multiple entries ##############
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

sort -k1,1 -k4,4 $resources/gencode.v43.primary_assembly.gene_name.bed | \
bedtools groupby -i stdin -g 1,4 -c 2,3 -o min,max | \
    awk -v OFS="\t" '{print $1,$3,$4,$2}' | \
    sort -k1,1 -k2,2n > \
    $resources/tmp && mv $resources/tmp $resources/gencode.v43.primary_assembly.gene_name.bed

bgzip $resources/gencode.v43.primary_assembly.gene_name.bed
tabix -p bed $resources/gencode.v43.primary_assembly.gene_name.bed.gz

###### Fix dbSNP contig names (1 -> chr1, etc.) #########
# Probably not best practice to mess around with VCF files but I'm not sure why the contigs were named like this,
# And I can't find a resource that has the correct naming scheme for hg38.
zcat $resources/00-common_all.vcf.gz | sed '/^[0-9XY]/s/^/chr/' | bgzip > $resources/dbsnp151_common.hg38.vcf.gz
tabix -p vcf $resources/dbsnp151_common.hg38.vcf.gz