rule bam_to_fastq:
    input:
        get_rna_sample_bams,
    output:
        fastq1="../results/{sample}/FASTQs/{sample}.R1.fastq",
        fastq2="../results/{sample}/FASTQs/{sample}.R2.fastq",
        singleton="../results/{sample}/FASTQs/{sample}.singleton.fastq",
        ambiguous="../results/{sample}/FASTQs/{sample}.ambiguous.fastq",
    log:
        "logs/picard/sam_to_fastq/{sample}.log",
    threads: 20
    shell:
        "samtools sort -n --threads 10 {input} | samtools fastq --threads 10 -1 {output.fastq1} -2 {output.fastq2} -s {output.singleton} -0 {output.ambiguous} -"