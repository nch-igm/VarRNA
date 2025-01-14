#!/bin/bash

resources="${1:-resources}"

# Reference data
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.p13.genome.fa.gz -P $resources
gunzip -c $resources/GRCh38.p13.genome.fa.gz > $resources/GRCh38.p13.genome.fa
gatk CreateSequenceDictionary -R $resources/GRCh38.p13.genome.fa
samtools faidx $resources/GRCh38.p13.genome.fa

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz -P $resources
gunzip -c $resources/gencode.v43.primary_assembly.annotation.gtf.gz > $resources/gencode.v43.primary_assembly.annotation.gtf

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.transcripts.fa.gz -P $resources

# Salmon
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz -P dependencies/
tar -xvzf dependencies/salmon-1.9.0_linux_x86_64.tar.gz -C dependencies/
dependencies/salmon-1.9.0_linux_x86_64/bin/salmon index -t $resources/gencode.v47.transcripts.fa.gz -i dependencies/gencode.v47.transcripts_index

# dbSNP known sites for GATK BaseRecalibrator
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz -P $resources

# Repeat Masker
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz -P $resources

#REDI Portal
wget http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz -P $resources

# ANNOVAR databases
dependencies/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar gnomad30_genome dependencies/human_db
dependencies/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar clinvar_20221231 dependencies/human_db
dependencies/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar dbnsfp42a dependencies/human_db
dependencies/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene dependencies/human_db
dependencies/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar cosmic70 dependencies/human_db



 
