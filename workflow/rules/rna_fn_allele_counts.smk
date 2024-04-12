rule create_fn_bed:
    input:
        rna="../results/{sample}/Matrix/{sample}.features.tsv",
        dna=expand("../results/{{sample}}/Matrix/{{sample}}_DNA_{variant_caller}.features.tsv", variant_caller=["HC", "MT2"]),
    output:
        "../results/{sample}/False_negative_allele_counts/{sample}.FalseNegative.bed",
    params:
        indir="../results/{sample}/Matrix",
    threads: 1
    shell:
        "python scripts/CreateFNAlleleCounts.py {wildcards.sample} {params.indir} {output}"

rule collect_allelic_counts:
    input:
        ref=config['ref'],
        bam="../results/{sample}/BAMs/{sample}.recal.bam",
        bed="../results/{sample}/False_negative_allele_counts/{sample}.FalseNegative.bed"
    output:
        "../results/{sample}/False_negative_allele_counts/{sample}.filtered.FalseNegative.AllelicCounts.tsv",
    threads: 1
    shell:
        "gatk CollectAllelicCounts "
        "-I {input.bam} "
        "-R {input.ref} "
        "-L {input.bed} "
        "-O {output}"