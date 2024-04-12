
sample=$1
input=$2
output_dir=$3
humandb=$4

../dependencies/annovar/table_annovar.pl --vcfinput $input $humandb --buildver hg38 --outfile $output_dir/$sample --remove \
    --protocol refGene,clinvar_20221231,gnomad30_genome,dbnsfp42a,cosmic70 \
    --operation g,f,f,f,f \
    --nastring . \
    --thread 10

mv $output_dir/$sample.hg38_multianno.vcf $output_dir/$sample.annovar.vcf
bgzip $output_dir/$sample.annovar.vcf
bcftools index --tbi --threads 20 $output_dir/$sample.annovar.vcf.gz