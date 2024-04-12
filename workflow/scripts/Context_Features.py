import sys, os
import pandas as pd
import itertools
from pyfaidx import Fasta


def linguistic_complexity(sequence):
    ''' 
    Get cartesian product for di and trinucleotide sequences.
    LC approximation described in this paper: https://doi.org/10.1016/S0097-8485(99)00007-8 
    '''
    U2_comb = itertools.product("ATCG", repeat=2)
    U3_comb = itertools.product("ATCG", repeat=3)
    U1 = 0
    for s in "ATCG":
        if s in sequence:
            U1 = U1 + 1

    U1 = min(1, U1/min(4, len(sequence))) # 4 for total number of nucleotides
    U2 = 0
    for s in U2_comb:
        s = s[0]+s[1]
        if s in sequence:
            U2 = U2 + 1

    U2 = min(1, U2/min(16, len(sequence)-1)) # 16 for total number of dinucleotide combinations
    U3 = 0
    for s in U3_comb:
        s = s[0]+s[1]+s[2]
        if s in sequence:
            U3 = U3 + 1
            
    U3 = min(1, U3/min(64, len(sequence)-2)) # ' '
    C = U1*U2*U3
    return(C)

def add_context_features(df, ref_genome):
    '''
    Add features based on sequencing context to dataframe of variants
    '''
    refdict=Fasta(ref_genome, as_raw=True)
    df['NUCLEOTIDE_CONTEXT'] = df['REF'] + df['ALT']

    nc_convert = {'AG':'transition', 'GA':'transition', 'CT':'transition', 'TC':'transition', 
    'AT':'transversion', 'TA':'transversion', 'AC':'transversion', 'CA':'transversion',
    'CG':'transversion', 'GC':'transversion', 'TG':'transversion', 'GT':'transversion'}
    df['NUCLEOTIDE_CONTEXT'] = df['NUCLEOTIDE_CONTEXT'].map(nc_convert)

    df['keytuple'] = list(zip(df['CHROM'], df['POS'].astype(int)))
    window_radius = 6
    df['locseq'] = df['keytuple'].apply(lambda x: refdict[x[0]][(x[1]-1-window_radius):x[1]+window_radius])
    df['locseq'] = df['locseq'].str[:window_radius] + df['locseq'].str[window_radius+1:]
    df['GC_CONTENT'] = (df["locseq"].str.count("G") + df["locseq"].str.count("C"))/(2*window_radius)
    df['LINGUISTIC_COMPLEXITY'] = df['locseq'].map(lambda x: linguistic_complexity(x))

    df = df.drop(['keytuple', 'locseq'], axis=1)
    return(df)