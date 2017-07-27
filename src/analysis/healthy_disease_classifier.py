#!/usr/bin/env python
"""
This script builds a general Healthy vs. Disease classifier combining
all datasets together. It builds classifiers with leave-one-disease out
and leave-one-dataset out cross validation.
"""
import argparse
import os
import sys
import pandas as pd
import numpy as np

from sklearn.metrics import roc_curve, auc

# Add util to the path
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
import FileIO as fio
from util import collapse_taxonomic_contents_df, prep_classifier, cv_and_roc

def test_dataset(rf, X_test, Y_test):
    """
    Test X_test with rf (which is a trained RandomForestClassifier).
    Returns fpr, tpr lists and roc_auc float.
    """

    # Test
    probs = rf.predict_proba(X_test)[:,1]
    predictions = rf.predict(X_test)

    fpr, tpr, thresholds = roc_curve(Y_test, probs)
    roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc


p = argparse.ArgumentParser()
p.add_argument('data_dir', help='path to directory with clean OTU tables and '
    'metadata.')
p.add_argument('out_file', help='file to write RF results to')
p.add_argument('--random-state', help='random state seed for classification',
    default=12345, type=int)
p.add_argument('--n-cv', help='number of cross validation folds [default: '
    + '%(default)s]', default=100, type=int)
args = p.parse_args()

datadir = args.data_dir
# Read in dfdict
dfdict = fio.read_dfdict_data(datadir)

## Collapse to genus level and relabel samples
for dataset in dfdict:
    # Collapse to genus level and relabel samples with dataset ID
    df = dfdict[dataset]['df']
    df = collapse_taxonomic_contents_df(df, 'genus')
    if dataset == 'edd_singh':
        df.index = ['cdi_singh-' + i for i in df.index]
    elif dataset == 'noncdi_schubert':
        df.index = ['cdi_schubert2-' + i for i in df.index]
    else:
        df.index = [dataset + '-' + i for i in df.index]
    dfdict[dataset]['df'] = df

    # Also relabel indices in metadata
    meta = dfdict[dataset]['meta']
    if dataset == 'edd_singh':
        meta.index = ['cdi_singh-' + i for i in meta.index]
    elif dataset == 'noncdi_schubert':
        meta.index = ['cdi_schubert2-' + i for i in meta.index]
    else:
        meta.index = [dataset + '-' + i for i in meta.index]
    dfdict[dataset]['meta'] = meta

## Concatenate OTU tables and corresponding metadata
# Only keep datasets with *healthy* controls
ignore_datasets = [d for d in dfdict
    if 'H' not in dfdict[d]['meta']['DiseaseState'].unique()]
# cdi_vincent only has five controls
ignore_datasets.append('cdi_vincent')

bigdf = pd.concat([dfdict[d]['df'] for d in dfdict if d not in ignore_datasets])
# Fill NaN's with zeros (i.e. unobserved OTUs)
bigdf = bigdf.fillna(0.0)
bigmeta = pd.concat([dfdict[d]['meta'] for d in dfdict
    if d not in ignore_datasets])

# excludes: 'postFMT_CDI', None, ' '
diseases = ['ASD', 'CD', 'CDI', 'nonCDI', 'CIRR', 'CRC', 'EDD', 'HIV',
            'MHE', 'NASH', 'OB',
            'PAR', 'PSA', 'RA', 'T1D', 'T2D', 'UC']
classes_list = [['H'], diseases]
[h_smpls, dis_smpls] = fio.get_samples(bigmeta, classes_list)

# Get rid of the ob_zhu healthy samples so they are not duplicated in the
# training set
h_smpls = [i for i in h_smpls if not i.startswith('ob_zhu')]

random_state = args.random_state

## Leave-one-dataset-out
datasets = list(set([i.split('-')[0] for i in bigdf.index]))
all_results = []
for d in datasets:
    print(d)
    # Train on all samples not in that dataset
    train_h = [i for i in h_smpls if not i.startswith(d)]
    # Get rid of the healthy samples that match with ob_zhu dataset
    if d == 'ob_zhu':
        train_h = [i for i in train_h if not i.startswith('nash_zhu')]

    train_dis = [i for i in dis_smpls if not i.startswith(d)]
    rf, X_train, Y_train = prep_classifier(
        bigdf, train_h, train_dis, random_state)
    rf = rf.fit(X_train, Y_train)

    # Test on that dataset
    test_h = [i for i in h_smpls if i.startswith(d)]
    # test_h will be empty for ob_zhu, the healthy samples are in nash_zhu
    if d == 'ob_zhu':
        test_h = [i for i in h_smpls if i.startswith('nash_zhu')]
    test_dis = [i for i in dis_smpls if i.startswith(d)]

    _, X_test, Y_test = prep_classifier(bigdf, test_h, test_dis, random_state)
    fpr, tpr, roc_auc = test_dataset(rf, X_test, Y_test)

    dis = d.split('_')[0] # note: edd_singh samples were relabeld to cdi_singh
    results = pd.DataFrame.from_dict(
        dict(zip(['dataset', 'disease', 'fpr', 'tpr', 'auc', 'classifier'],
                 (d, dis, fpr, tpr, roc_auc, 'dataset_out'))))
    all_results.append(results)

## Leave-one-disease-out
diseases = list(set([i.split('_')[0] for i in bigdf.index]))
for d in diseases:
    print(d)
    train_h = [i for i in h_smpls if not i.startswith(d)]

    # For ob comparison, get rid of the healthy controls from nash_zhu
    # because they belong to ob_zhu dataset too
    if d == 'ob':
        train_h = [i for i in train_h if not i.startswith('nash_zhu')]

    train_dis = [i for i in dis_smpls if not i.startswith(d)]
    rf, X_train, Y_train = prep_classifier(
        bigdf, train_h, train_dis, random_state)

    test_h = [i for i in h_smpls if i.startswith(d)]
    # Add the ob_zhu healthy samples here
    if d == 'ob':
        test_h += [i for i in h_smpls if i.startswith('nash_zhu')]

    test_dis = [i for i in dis_smpls if i.startswith(d)]
    _, X_test, Y_test = prep_classifier(bigdf, test_h, test_dis, random_state)

    # Train
    rf = rf.fit(X_train, Y_train)
    # Test
    probs = rf.predict_proba(X_test)[:,1]

    # Get the fpr, tpr, and AUC for disease-wise results
    fpr, tpr, thresholds = roc_curve(Y_test, probs)
    roc_auc = auc(fpr, tpr)

    results = pd.DataFrame.from_dict(
        dict(zip(['disease', 'fpr', 'tpr', 'auc', 'classifier'],
                 (d, fpr, tpr, roc_auc, 'disease_out'))))
    all_results.append(results)

    # Also get the dataset-wise fpr, tpr, and AUC
    # Make dataframe with probabilities and Y_test
    all_smpls = test_h + test_dis
    probs_df = pd.DataFrame(probs, index=all_smpls, columns=['probs'])
    probs_df['y_test'] = Y_test
    probs_df['dataset'] = [i.split('-')[0] for i in probs_df.index]
    if d == 'ob':
        probs_df = probs_df.replace('nash_zhu', 'ob_zhu')

    for g, subdf in probs_df.groupby('dataset'):
        fpr, tpr, thresholds = roc_curve(subdf['y_test'], subdf['probs'])
        roc_auc = auc(fpr, tpr)

        results = pd.DataFrame.from_dict(
            dict(zip(['dataset', 'disease', 'fpr', 'tpr', 'auc', 'classifier'],
                     (g, d, fpr, tpr, roc_auc, 'disease_out'))))
        all_results.append(results)

all_results_df = pd.concat(all_results)
all_results_df.to_csv(args.out_file, sep='\t', index=False)
