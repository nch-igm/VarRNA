import sys, os
import pandas as pd

genes_bed = sys.argv[1]
salmon = sys.argv[2]
output=sys.argv[3]


genes = pd.read_csv(genes_bed, compression='gzip', header=None, sep='\t')

salmon = pd.read_csv(salmon, sep='\t')
salmon['GeneName'] = salmon['Name'].str.split('|').str[5]
salmon['Type'] = salmon['Name'].str.split('|').str[7]
salmon = salmon.sort_values('TPM', ascending=False).groupby(['GeneName'], sort=False).first().reset_index()

genes = pd.read_csv(genes_bed, compression='gzip', header=None, sep='\t')
genes.columns = ['chrom', 'start', 'stop', 'GeneName']

df = genes.merge(salmon[['GeneName', 'TPM', 'Length', 'EffectiveLength', 'NumReads', 'Type']], how='left')
df.to_csv(output, index=False, header=None, sep='\t')


