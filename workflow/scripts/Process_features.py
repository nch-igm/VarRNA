import sys
import pandas as pd
import numpy as np
from pyfaidx import Fasta
from  Context_Features import *

input = sys.argv[1]
output=sys.argv[2]
ref_genome=sys.argv[3]

# List of features for machine learning in order
feature_columns = ['QUAL', 'AF', 'GQ', 'INFO/BaseQRankSum', 'INFO/FS', 'INFO/MPOS',
       'INFO/QD', 'INFO/ReadPosRankSum', 'INFO/SOR', 'INFO/AF_raw',
       'INFO/LENGTH', 'AD_REF', 'AD_ALT', 'MBQ_REF', 'MBQ_ALT', 'MFRL_REF',
       'MFRL_ALT', 'GC_CONTENT', 'LINGUISTIC_COMPLEXITY',
       'INFO/TRANS_TYPE_non_stop_decay',
       'INFO/TRANS_TYPE_nonsense_mediated_decay',
       'INFO/TRANS_TYPE_processed_transcript',
       'INFO/TRANS_TYPE_protein_coding', 'INFO/TRANS_TYPE_retained_intron',
       'NUCLEOTIDE_CONTEXT_transition', 'NUCLEOTIDE_CONTEXT_transversion',
       'INFO/SIFT_pred_D', 'INFO/SIFT_pred_T', 'INFO/SIFT4G_pred_D',
       'INFO/SIFT4G_pred_T', 'INFO/Polyphen2_HDIV_pred_B',
       'INFO/Polyphen2_HDIV_pred_D', 'INFO/Polyphen2_HDIV_pred_P',
       'INFO/Polyphen2_HVAR_pred_B', 'INFO/Polyphen2_HVAR_pred_D',
       'INFO/Polyphen2_HVAR_pred_P', 'INFO/LRT_pred_D', 'INFO/LRT_pred_N',
       'INFO/LRT_pred_U', 'INFO/MutationTaster_pred_A',
       'INFO/MutationTaster_pred_D', 'INFO/MutationTaster_pred_N',
       'INFO/MutationTaster_pred_P', 'INFO/FATHMM_pred_D',
       'INFO/FATHMM_pred_T', 'INFO/PROVEAN_pred_D', 'INFO/PROVEAN_pred_N',
       'INFO/MetaSVM_pred_D', 'INFO/MetaLR_pred_D', 'INFO/MetaRNN_pred_D',
       'INFO/PrimateAI_pred_D', 'INFO/PrimateAI_pred_T', 'INFO/DEOGEN2_pred_D',
       'INFO/DEOGEN2_pred_T', 'INFO/BayesDel_addAF_pred_D',
       'INFO/BayesDel_noAF_pred_D', 'INFO/ClinPred_pred_D',
       'INFO/ClinPred_pred_T', 'INFO/phyloP100way_vertebrate',
       'INFO/phyloP100way_vertebrate_rankscore', 'INFO/phyloP30way_mammalian',
       'INFO/phyloP30way_mammalian_rankscore', 'INFO/cosmic70', 'indel']


df = pd.read_csv(input, sep='\t')

df.columns = df.columns.str.lstrip("%")
df = df.replace('.', np.NaN)
df = df.apply(pd.to_numeric, errors='ignore')
df['POS'] = df['POS'].astype(str)

# Separate comma separated values
df[['AD_REF','AD_ALT']] = df['AD'].str.split(',',expand=True).astype(int)
df[['MBQ_REF', 'MBQ_ALT']] = df['INFO/MBQ'].str.split(',',expand=True).astype(int)
df[['MFRL_REF', 'MFRL_ALT']] = df['INFO/MFRL'].str.split(',',expand=True).astype(int)
df = df.drop(['AD', 'INFO/MBQ', 'INFO/MFRL'], axis=1)

# Convert gnomad columns missing values to zero
df['INFO/AF_raw'] = df['INFO/AF_raw'].replace(np.NAN, 0)

# Add context features
nc_convert = {'AG':'transition', 'GA':'transition', 'CT':'transition', 'TC':'transition', 
    'AT':'transversion', 'TA':'transversion', 'AC':'transversion', 'CA':'transversion',
    'CG':'transversion', 'GC':'transversion', 'TG':'transversion', 'GT':'transversion'}

df['NUCLEOTIDE_CONTEXT'] = df['REF'] + df['ALT']
df['NUCLEOTIDE_CONTEXT'] = df['NUCLEOTIDE_CONTEXT'].map(nc_convert)

refdict=Fasta(ref_genome, as_raw=True)

df['keytuple'] = list(zip(df['CHROM'], df['POS'].astype(int)))
window_radius = 6
df['locseq'] = df['keytuple'].apply(lambda x: refdict[x[0]][(x[1]-1-window_radius):x[1]+window_radius])
df['locseq'] = df['locseq'].str[:window_radius] + df['locseq'].str[window_radius+1:]

df['GC_CONTENT'] = (df["locseq"].str.count("G") + df["locseq"].str.count("C"))/(2*window_radius)
df['LINGUISTIC_COMPLEXITY'] = df['locseq'].map(lambda x: linguistic_complexity(x))

df['indel'] = np.where((df['REF'].str.len() > 1) | (df['ALT'].str.len() > 1 ), 1, 0)

df = df.drop(['keytuple', 'locseq'], axis=1)

# Get dummy values for the categorical columns
df = pd.get_dummies(df, columns=['INFO/ExonicFunc.refGene', 'INFO/CLNREVSTAT', 'INFO/CLNSIG', 'INFO/TRANS_TYPE', 'NUCLEOTIDE_CONTEXT',
                        'INFO/SIFT_pred', 'INFO/SIFT4G_pred', 'INFO/Polyphen2_HDIV_pred',
                        'INFO/Polyphen2_HVAR_pred', 'INFO/LRT_pred', 'INFO/MutationTaster_pred',
                        'INFO/MutationAssessor_pred', 'INFO/FATHMM_pred', 'INFO/PROVEAN_pred',
                        'INFO/MetaSVM_pred', 'INFO/MetaLR_pred', 'INFO/MetaRNN_pred',
                        'INFO/PrimateAI_pred', 'INFO/DEOGEN2_pred', 'INFO/BayesDel_addAF_pred',
                        'INFO/BayesDel_noAF_pred', 'INFO/ClinPred_pred'])


# Row filtering
df = df[(df['DP'] >=10)]
df = df[df['QUAL'] >= 100]
df = df[df['INFO/AF_raw'] <= .001]

# Subset columns to match features for machine learning
df.index = df['CHROM'] + '_' + df['POS'].astype(str) + '_' + df['REF'] + '_' + df['ALT'] + '_' + df['SAMPLE']

cols_to_add = list(set(feature_columns) - set(df.columns))
df[cols_to_add] = 0
df = df[feature_columns]

df.to_csv(output)

