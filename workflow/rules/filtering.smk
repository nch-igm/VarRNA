rule filter_vcfs_bed_rna:
    input:
        rep_mask = config["annotation_resources"]["rep_mask"],
        exome_bed = config["annotation_resources"]["exome_bed"],
        vcf = "../results/{sample}/VCFs/{sample}.vcf.gz",
        vcf_tbi="../results/{sample}/VCFs/{sample}.vcf.gz.tbi",
    output:
        vcf="../results/{sample}/VCFs/filtered/{sample}.filtered.vcf.gz",
        vcf_tbi="../results/{sample}/VCFs/filtered/{sample}.filtered.vcf.gz.tbi",
    threads: 1
    shell:
        "bash scripts/FilterVCF.sh {input.vcf} {output.vcf} {input.rep_mask} {input.exome_bed}"
