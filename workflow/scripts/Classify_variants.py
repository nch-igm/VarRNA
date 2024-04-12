#!/usr/bin/env python3

import sys, os
import argparse
import pickle
import pandas as pd

input=sys.argv[1]
output=sys.argv[2]
input_tsv=sys.argv[3]
model_ta=sys.argv[4]
model_gs=sys.argv[5]

def predict_labels(input_df, model_path, label_dict):
    with open(model_path, 'rb') as f:
        xgb_search = pickle.load(f)
    clf = xgb_search.best_estimator_

    y_hat = clf.predict(input_df)
    y_hat = list(map(label_dict.get, y_hat))
    return(y_hat)


def combine_labels(y_hat_ta, y_hat_gs):
    y_hat_comb = []
    for i,label in enumerate(y_hat_ta):
        if label == 'True Variant':
            val = y_hat_gs[i]
            y_hat_comb.append(val)
        else:
            y_hat_comb.append(label)
    return(y_hat_comb)


df = pd.read_csv(input, index_col=0)

# Predict variant labels
label_dicts = [{1: 'True Variant', 0: 'Artifact'}, {0: 'Germline', 1: 'Somatic'}]
y_hat_ta = predict_labels(df, model_ta, label_dicts[0])
y_hat_gs = predict_labels(df, model_gs, label_dicts[1])
y_hat = combine_labels(y_hat_ta, y_hat_gs)

# Add label to annotated variants tsv
df_out = pd.read_csv(input_tsv, sep='\t')
df_out.columns = df_out.columns.str.lstrip("%")

df_out.index = df_out['CHROM'] + '_' + df_out['POS'].astype(str) + '_' + df_out['REF'] + '_' + df_out['ALT'] + '_' + df_out['SAMPLE']
df_out = df_out.loc[df.index]

df_out['Model_Label'] = y_hat

df_out.to_csv(output)

