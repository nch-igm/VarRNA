rule classify_variants:
    input:
        feats_csv="../results/{sample}/Matrix/{sample}.processed_features.csv",
        ann_tsv="../results/{sample}/Matrix/{sample}.features.tsv",
        model_ta="data/xgb_trueartifact.pkl",
        model_gs="data/xgb_germsom.pkl",
    output:
        "../results/{sample}/Predictions/{sample}.annotated_predictions.csv"
    threads: 1
    shell:
        "python scripts/Classify_variants.py {input.feats_csv} {output} {input.ann_tsv} {input.model_ta} {input.model_gs}"