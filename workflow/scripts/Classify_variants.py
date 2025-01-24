#!/usr/bin/env python3

import sys
import pickle
import pandas as pd
import numpy as np

input=sys.argv[1]
output=sys.argv[2]
input_tsv=sys.argv[3]
model_ta=sys.argv[4]
model_gs=sys.argv[5]

def predict_probs(input_df, model_path):
    """
    Predict the probabilities of the input data using a pre-trained model.

    Parameters:
    input_df (pandas.DataFrame): The input data for which probabilities need to be predicted.
    model_path (str): The file path to the pre-trained model (in pickle format).

    Returns:
    numpy.ndarray: An array of predicted probabilities for each class.
    """
    with open(model_path, 'rb') as f:
        xgb_search = pickle.load(f)
    clf = xgb_search.best_estimator_

    y_prob = clf.predict_proba(input_df)
    return(y_prob)


def combine_labels(y_prob_ta, y_prob_gs, label_dict):
    """
    Combines TA (True/Artifact) and GS (Germline/Somatic) labels based on their probabilities.

    Args:
        y_prob_ta (list of list of float): Probabilities for TA classification.
        y_prob_gs (list of list of float): Probabilities for GS classification.
        label_dict (dict): Dictionary mapping indices to labels for TA and GS.

    Returns:
        tuple: A tuple containing:
            - y_hat_comb (list of str): Combined labels.
            - y_prob_comb (list of float): Combined probabilities.
    """
    y_hat_comb = []
    y_prob_comb = []

    # Map TA labels and GS labels
    y_hat_ta = [label_dict[0][np.argmax(probs)] for probs in y_prob_ta]

    for ta_label, ta_prob, gs_prob in zip(y_hat_ta, y_prob_ta, y_prob_gs):
        if ta_label == 'True Variant':
            # Scale probabilities by the True Variant probability
            scaled_germline_prob = gs_prob[0] * ta_prob[1]
            scaled_somatic_prob = gs_prob[1] * ta_prob[1]

            # Choose Germline or Somatic based on scaled probabilities
            if scaled_germline_prob > scaled_somatic_prob:
                y_hat_comb.append('Germline')
                y_prob_comb.append(scaled_germline_prob)
            else:
                y_hat_comb.append('Somatic')
                y_prob_comb.append(scaled_somatic_prob)
        else:
            # For Artifact, use existing probability
            y_hat_comb.append(ta_label)
            y_prob_comb.append(ta_prob[0])

    return(y_hat_comb, y_prob_comb)


df = pd.read_csv(input, index_col=0)

# Predict variant labels
label_dicts = [{1: 'True Variant', 0: 'Artifact'}, {0: 'Germline', 1: 'Somatic'}]
y_prob_ta = predict_probs(df, model_ta)
df_gs = df.drop(columns=['AF'])
y_prob_gs = predict_probs(df_gs, model_gs)
y_hat_comb, y_prob_comb = combine_labels(y_prob_ta, y_prob_gs, label_dicts)

# Add label and probability to annotated variants tsv
df_out = pd.read_csv(input_tsv, sep='\t')
df_out.columns = df_out.columns.str.lstrip("%")

df_out.index = df_out['CHROM'] + '_' + df_out['POS'].astype(str) + '_' + df_out['REF'] + '_' + df_out['ALT'] + '_' + df_out['SAMPLE']
df_out = df_out.loc[df.index]

df_out['Model_Label'] = y_hat_comb
df_out['Model_Prob'] = y_prob_comb

df_out.to_csv(output)