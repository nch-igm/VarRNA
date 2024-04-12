rule process_features:
    input:
        tsv="../results/{sample}/Matrix/{sample}.features.tsv",
        ref=config["reference"]["fasta"],
    output:
        "../results/{sample}/Matrix/{sample}.processed_features.csv",
    threads: 1
    shell:
        "python scripts/Process_features.py {input.tsv} {output} {input.ref}"