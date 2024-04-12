import sys, os
import pandas as pd

sample = sys.argv[1]
indir = sys.argv[2]
fout = sys.argv[3]

columns = ['%CHROM', '%POS', '%REF', '%ALT', '%SAMPLE', '%QUAL', '%GT']

rna_results = pd.read_csv(f'{indir}/{sample}.features.tsv', sep='\t', usecols=columns)
rna_results['RNA_variant_called'] = 'T'

hc_results = pd.read_csv(f'{indir}/{sample}_DNA_HC.features.tsv', sep='\t', usecols=columns)
mt2_results = pd.read_csv(f'{indir}/{sample}_DNA_MT2.features.tsv', sep='\t', usecols=columns)
dna_results = pd.concat([hc_results, mt2_results])

df = rna_results.merge(dna_results, on=['%CHROM', '%POS', '%REF', '%ALT'], how='outer')

df_fn = df[df['RNA_variant_called'] != 'T']

out_bed = df_fn[['%CHROM', '%POS']]
out_bed = out_bed.drop_duplicates()
out_bed['POS1'] = out_bed['%POS']-1
out_bed[['%CHROM', 'POS1', '%POS']].to_csv(fout, header=False, index=False, sep='\t')