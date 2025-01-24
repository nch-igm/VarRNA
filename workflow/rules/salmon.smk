rule salmon_quant_reads:
    input:
        r1="../results/{sample}/FASTQs/{sample}.R1.fastq.gz",
        r2="../results/{sample}/FASTQs/{sample}.R2.fastq.gz",
        index=config["dependencies"]["salmon_index"],
    output:
        quant="../results/{sample}/salmon/quant.sf",
        lib="../results/{sample}/salmon/lib_format_counts.json",
    params:
        salmon=config["dependencies"]["salmon"]
    log:
        "logs/salmon/{sample}.log",
    threads: 10
    shell:
        "{params.salmon} quant --threads 10 -l A -o ../results/{wildcards.sample}/salmon/ -i {input.index} -1 {input.r1} -2 {input.r2}"

# rule modify_names:
#    input:
#        quant="../results/{sample}/salmon/quant.sf",
#        gentrome=config["reference"]["gentrome"],
#    output:
#        "../results/{sample}/salmon/quant.sf.orig",
#    log:
#        "logs/salmon/modify_names.{sample}.log",
#    shell:
#        "scripts/map_ids.pl {input.quant} {input.gentrome}"