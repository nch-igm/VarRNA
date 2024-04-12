import sys, os
import pandas as pd
import numpy as np
import itertools
from pyfaidx import Fasta
from  Context_Features import *

input = sys.argv[1]
output=sys.argv[2]
ref_genome=sys.argv[3]

# List of features for machine learning in order
feature_columns = ['QUAL', 'AF', 'DP', 'GQ', 'INFO/AC', 'INFO/AF', 'INFO/BaseQRankSum', 'INFO/DP', 'INFO/FS', 'INFO/MLEAC', 'INFO/MLEAF', 'INFO/MPOS',
       'INFO/NCount', 'INFO/QD', 'INFO/ReadPosRankSum', 'INFO/SOR', 'INFO/AF.1', 'INFO/AF_raw', 'INFO/AF_male', 'INFO/AF_female',
       'INFO/AF_afr', 'INFO/AF_ami', 'INFO/AF_amr', 'INFO/AF_asj', 'INFO/AF_eas', 'INFO/AF_fin', 'INFO/AF_nfe', 'INFO/AF_oth',
       'INFO/AF_sas', 'INFO/TPM', 'INFO/LENGTH', 'INFO/EFF_LENGTH', 'INFO/NUM_READS', 'AD_REF', 'AD_ALT', 'MBQ_REF', 'MBQ_ALT', 'MFRL_REF',
       'MFRL_ALT', 'GC_CONTENT', 'LINGUISTIC_COMPLEXITY', 'INFO/ExonicFunc.refGene_frameshift_deletion', 'INFO/ExonicFunc.refGene_frameshift_insertion',
       'INFO/ExonicFunc.refGene_nonframeshift_deletion', 'INFO/ExonicFunc.refGene_nonframeshift_insertion', 'INFO/ExonicFunc.refGene_nonsynonymous_SNV',
       'INFO/ExonicFunc.refGene_startloss', 'INFO/ExonicFunc.refGene_stopgain', 'INFO/ExonicFunc.refGene_stoploss', 'INFO/ExonicFunc.refGene_synonymous_SNV',
       'INFO/ExonicFunc.refGene_unknown', 'INFO/CLNREVSTAT_criteria_provided,_conflicting_interpretations',
       'INFO/CLNREVSTAT_criteria_provided,_multiple_submitters,_no_conflicts', 'INFO/CLNREVSTAT_criteria_provided,_single_submitter',
       'INFO/CLNREVSTAT_no_assertion_criteria_provided', 'INFO/CLNREVSTAT_no_assertion_provided', 'INFO/CLNREVSTAT_reviewed_by_expert_panel', 'INFO/CLNSIG_Benign',
       'INFO/CLNSIG_Benign/Likely_benign', 'INFO/CLNSIG_Conflicting_interpretations_of_pathogenicity',
       'INFO/CLNSIG_Conflicting_interpretations_of_pathogenicity|risk_factor', 'INFO/CLNSIG_Likely_benign', 'INFO/CLNSIG_Likely_pathogenic',
       'INFO/CLNSIG_Pathogenic', 'INFO/CLNSIG_Pathogenic/Likely_pathogenic', 'INFO/CLNSIG_Uncertain_significance', 'INFO/CLNSIG_association',
       'INFO/CLNSIG_not_provided', 'INFO/TRANS_TYPE_TEC', 'INFO/TRANS_TYPE_antisense', 'INFO/TRANS_TYPE_non_stop_decay',
       'INFO/TRANS_TYPE_nonsense_mediated_decay', 'INFO/TRANS_TYPE_processed_transcript',
       'INFO/TRANS_TYPE_protein_coding', 'INFO/TRANS_TYPE_retained_intron', 'NUCLEOTIDE_CONTEXT_transition', 'NUCLEOTIDE_CONTEXT_transversion',
       'INFO/SIFT_pred_D', 'INFO/SIFT_pred_T', 'INFO/SIFT4G_pred_D', 'INFO/SIFT4G_pred_T', 'INFO/Polyphen2_HDIV_pred_B', 'INFO/Polyphen2_HDIV_pred_D', 
       'INFO/Polyphen2_HDIV_pred_P', 'INFO/Polyphen2_HVAR_pred_B', 'INFO/Polyphen2_HVAR_pred_D',
       'INFO/Polyphen2_HVAR_pred_P', 'INFO/LRT_pred_D', 'INFO/LRT_pred_N', 'INFO/LRT_pred_U', 'INFO/MutationTaster_pred_A',
       'INFO/MutationTaster_pred_D', 'INFO/MutationTaster_pred_N', 'INFO/MutationTaster_pred_P', 'INFO/MutationAssessor_pred_H',
       'INFO/MutationAssessor_pred_L', 'INFO/MutationAssessor_pred_M', 'INFO/MutationAssessor_pred_N', 'INFO/FATHMM_pred_D',
       'INFO/FATHMM_pred_T', 'INFO/PROVEAN_pred_D', 'INFO/PROVEAN_pred_N', 'INFO/MetaSVM_pred_D', 'INFO/MetaSVM_pred_T', 'INFO/MetaLR_pred_D',
       'INFO/MetaLR_pred_T', 'INFO/MetaRNN_pred_D', 'INFO/MetaRNN_pred_T', 'INFO/PrimateAI_pred_D', 'INFO/PrimateAI_pred_T', 'INFO/DEOGEN2_pred_D',
       'INFO/DEOGEN2_pred_T', 'INFO/BayesDel_addAF_pred_D', 'INFO/BayesDel_addAF_pred_T', 'INFO/BayesDel_noAF_pred_D',
       'INFO/BayesDel_noAF_pred_T', 'INFO/ClinPred_pred_D', 'INFO/ClinPred_pred_T', 'INFO/phyloP100way_vertebrate', 'INFO/phyloP100way_vertebrate_rankscore',
       'INFO/phyloP30way_mammalian', 'INFO/phyloP30way_mammalian_rankscore', 'INFO/cosmic70', 'indel']


df = pd.read_csv(input, sep='\t')

df.columns = df.columns.str.lstrip("%")
df = df.replace('.', np.NaN)
df = df.apply(pd.to_numeric, errors='ignore')
df['POS'] = df['POS'].astype(str)

# Separate comma separated values
df[['AD_REF','AD_ALT']] = df['AD'].str.split(',',expand=True).astype(int)
df[['MBQ_REF', 'MBQ_ALT']] = df['INFO/MBQ'].str.split(',',expand=True).astype(int)
df[['MFRL_REF', 'MFRL_ALT']] = df['INFO/MFRL'].str.split(',',expand=True).astype(int)
df[['MMQ_REF', 'MMQ_ALT']] = df['INFO/MMQ'].str.split(',',expand=True).astype(int)
df = df.drop(['AD', 'INFO/MBQ', 'INFO/MFRL', 'INFO/MMQ'], axis=1)

# Convert gnomad columns missing values to zero
gnomad_cols = ['INFO/AF.1', 'INFO/AF_raw', 'INFO/AF_male', 'INFO/AF_female', 'INFO/AF_afr', 'INFO/AF_ami', 'INFO/AF_amr', 'INFO/AF_asj', 'INFO/AF_eas', 'INFO/AF_fin', 'INFO/AF_nfe', 'INFO/AF_oth', 'INFO/AF_sas']
df[gnomad_cols] = df[gnomad_cols].replace(np.NAN, 0)

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
                        'INFO/BayesDel_noAF_pred', 'INFO/ClinPred_pred', ])

# Drop columns with zero variance (e.g. all values are the same)
# Also drop GT since it is exactly the same info as INFO/AC
df = df.drop(['INFO/AN', 'INFO/ExcessHet', 'INFO/MQ', 'INFO/MQRankSum', 'MMQ_REF', 'MMQ_ALT', 'GT'], axis=1)

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

