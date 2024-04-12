rule query_vcf:
    input:
        vcf = "../results/{sample}/VCFs/annotated/{sample}.annotated.vcf.gz",
        vcf_tbi = "../results/{sample}/VCFs/annotated/{sample}.annotated.vcf.gz.tbi",
    output:
        "../results/{sample}/Matrix/{sample}.features.tsv"
    params:
        querying = get_bcftools_query_values(),
    threads: 1
    shell:
        """
        echo -e '{params.querying}' > {output}
        bcftools query -f '[{params.querying}\n]' {input.vcf} >> {output}
        """