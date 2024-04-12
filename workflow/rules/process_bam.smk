rule add_read_groups:
    input:
        get_rna_sample_bams,
    output:
        temp("../results/{sample}/BAMs/{sample}.AddReadGroup.bam"),
    log:
        "logs/picard/addreadgroup/{sample}.log",
    threads: 6
    shell:
        "gatk AddOrReplaceReadGroups "
         "--I {input} "
         "--O {output} "
         "--RGID 1 "
         "--RGLB lib1 "
         "--RGPL ILLUMINA "
         "--RGPU unit1 "
         "--RGSM {wildcards.sample} "
         "--SORT_ORDER coordinate"

rule mark_duplicates:
    input:
        bams="../results/{sample}/BAMs/{sample}.AddReadGroup.bam",
    output:
        bam=temp("../results/{sample}/BAMs/{sample}.marked_duplicates.bam"),
        metrics="../results/{sample}/BAMs/{sample}.metrics.txt",
    log:
        "logs/picard/Mark_duplicates/{sample}.log",
    threads: 8
    shell:
         "gatk MarkDuplicates "
         "-I {input} "
         "-O {output.bam} "
         "-M {output.metrics} "
         "--CREATE_INDEX true"

rule splitncigarreads:
    input:
        bam="../results/{sample}/BAMs/{sample}.marked_duplicates.bam",
        ref=config["reference"]["fasta"],
    output:
        bam=temp("../results/{sample}/BAMs/{sample}.splitncigar.bam"),
        bai=temp("../results/{sample}/BAMs/{sample}.splitncigar.bai"),
    log:
        "logs/gatk/splitNCIGARreads/{sample}.log",
    threads: 8
    shell:
         "gatk SplitNCigarReads "
         "-R {input.ref} "
         "-I {input.bam} "
         "-O {output.bam}"

rule gatk_baserecalibrator:
    input:
        bam="../results/{sample}/BAMs/{sample}.splitncigar.bam",
        ref=config["reference"]["fasta"],
        dbsnp=config["annotation_resources"]["dbsnp"],
    output:
        recal_table="../results/{sample}/BAMs/{sample}.recal.grp",
    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    threads: 6
    shell:
         "gatk BaseRecalibrator "
         "-I {input.bam} "
         "-R {input.ref} "
         "--known-sites {input.dbsnp} "
         "-O {output}"

rule gatk_applybqsr:
    input:
        bam="../results/{sample}/BAMs/{sample}.splitncigar.bam",
        ref=config["reference"]["fasta"],
        recal_file="../results/{sample}/BAMs/{sample}.recal.grp",
    output:
        bam="../results/{sample}/BAMs/{sample}.recal.bam",
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log",
    threads: 6
    shell:
         "gatk ApplyBQSR "
         "-R {input.ref} "
         "-I {input.bam} "
         "--bqsr-recal-file {input.recal_file} "
         "-O {output} "
         "--create-output-bam-index"