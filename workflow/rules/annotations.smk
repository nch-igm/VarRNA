rule norm:
    input:
        vcf="../results/{sample}/VCFs/filtered/{sample}.filtered.vcf.gz",
        vcf_tbi="../results/{sample}/VCFs/filtered/{sample}.filtered.vcf.gz.tbi",
        ref=config["reference"]["fasta"],
    output:
        temp("../results/{sample}/VCFs/annotated/{sample}.norm.vcf.gz"),
    threads: 10
    shell:
        "bcftools norm -f {input.ref} -a -m - -Oz -o {output} {input.vcf}"

rule norm_tbi:
    input:
        "../results/{sample}/VCFs/annotated/{sample}.norm.vcf.gz",
    output:
        temp("../results/{sample}/VCFs/annotated/{sample}.norm.vcf.gz.tbi")
    threads: 10
    shell:
        "bcftools index --tbi --threads 20 {input}"

rule annovar:
    input:
        vcf="../results/{sample}/VCFs/annotated/{sample}.norm.vcf.gz",
        vcf_tbi="../results/{sample}/VCFs/annotated/{sample}.norm.vcf.gz.tbi",
    output:
        vcf=temp("../results/{sample}/VCFs/annotated/{sample}.annovar.vcf.gz"),
        vcf_tbi=temp("../results/{sample}/VCFs/annotated/{sample}.annovar.vcf.gz.tbi"),
    params:
        human_db=config["dependencies"]["annovar_db"],
        output_dir="../results/{sample}/VCFs/annotated/",
    threads: 10
    shell:
        "bash scripts/Annovar.sh {wildcards.sample} {input.vcf} {params.output_dir} {params.human_db}"

rule combine_gene_tpm:
    input:
        genes_bed=config["annotation_resources"]["genes_bed"],
        salmon="../results/{sample}/salmon/quant.sf",
        #salmon_mod="../results/{sample}/salmon/quant.sf.orig",
    output:
        temp("../results/{sample}/salmon/genes_tpm.bed"),
    threads: 1
    shell:
        "python scripts/Combine_genes_tpm.py {input.genes_bed} {input.salmon} {output}"

rule bgzip_gene_tpm:
    input:
        "../results/{sample}/salmon/genes_tpm.bed",
    output:
        "../results/{sample}/salmon/genes_tpm.bed.gz",
    threads: 1
    shell:
        "bgzip {input}"

rule gene_tpm_tbi:
    input:
        "../results/{sample}/salmon/genes_tpm.bed.gz",
    output:
        temp("../results/{sample}/salmon/genes_tpm.bed.gz.tbi"),
    threads: 1
    shell:
        "tabix -p bed {input}"

rule annot_genes:
    input:
        vcf="../results/{sample}/VCFs/annotated/{sample}.annovar.vcf.gz",
        vcf_tbi="../results/{sample}/VCFs/annotated/{sample}.annovar.vcf.gz.tbi",
        genes_tpm="../results/{sample}/salmon/genes_tpm.bed.gz",
        genes_tpm_tbi="../results/{sample}/salmon/genes_tpm.bed.gz.tbi",
        annots_header="scripts/annots.hdr",
    output:
        "../results/{sample}/VCFs/annotated/{sample}.annotated.vcf.gz",
    threads: 10
    shell:
        "bcftools annotate --threads 20 -a {input.genes_tpm} -h {input.annots_header} -c CHROM,FROM,TO,GENE,TPM,LENGTH,EFF_LENGTH,NUM_READS,TRANS_TYPE -Oz -o {output} {input.vcf}"

rule annot_genes_tbi:
    input:
        "../results/{sample}/VCFs/annotated/{sample}.annotated.vcf.gz",
    output:
        temp("../results/{sample}/VCFs/annotated/{sample}.annotated.vcf.gz.tbi"),
    threads: 10
    shell:
        "bcftools index --tbi --threads 20 {input}"



    
