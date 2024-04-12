rule haplotype_caller:
    input:
        bam="../results/{sample}/BAMs/{sample}.recal.bam",
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("../results/{sample}/VCFs/{sample}.{chrom}.vcf.gz"),
        vcf_tbi=temp("../results/{sample}/VCFs/{sample}.{chrom}.vcf.gz.tbi"),
    log:
        "logs/gatk/haplotypecaller/{sample}.{chrom}.log",
    threads: 4
    shell:
        "gatk HaplotypeCaller "
        "-I {input.bam} "
        "-R {input.ref} "
        "-O {output.vcf} "
        "-L {wildcards.chrom} "
        "--dont-use-soft-clipped-bases "
        "--standard-min-confidence-threshold-for-calling 20 "
        "--max-reads-per-alignment-start 0 "
        "-A BaseQuality "
        "-A BaseQualityRankSumTest "
        "-A FragmentLength "
        "-A MappingQualityRankSumTest "
        "-A MappingQuality "
        "-A AlleleFraction "
        "-A CountNs "
        "-A Coverage "
        "-A ReadPosition "
        "-A ReadPosRankSumTest"

rule concat_vcfs:
    input:
        expand("../results/{{sample}}/VCFs/{{sample}}.{chrom}.vcf.gz", chrom=get_chroms()),
    output:
        "../results/{sample}/VCFs/{sample}.vcf.gz",
        "../results/{sample}/VCFs/{sample}.vcf.gz.tbi",
    log:
        "logs/bcftools/{sample}.concat.log"
    params:
        vcf_dir = "../results/{sample}/VCFs",
    threads: 10
    shell:
        "bash scripts/ConcatVCFs.sh {wildcards.sample} {params.vcf_dir}"
