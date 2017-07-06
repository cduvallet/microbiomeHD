#!/usr/bin/env python
"""
This script builds a general Healthy vs. Disease classifier combining
all datasets together.
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
    df.index = [dataset + '-' + i for i in df.index]
    dfdict[dataset]['df'] = df

    # Also relabel indices in metadata
    meta = dfdict[dataset]['meta']
    meta.index = [dataset + '-' + i for i in meta.index]
    dfdict[dataset]['meta'] = meta

## Concatenate OTU tables and corresponding metadata
# Only keep datasets with *healthy* controls
ignore_datasets = [d for d in dfdict
    if 'H' not in dfdict[d]['meta']['DiseaseState'].unique()]
# This data is a duplicate of nash_zhu
ignore_datasets += ['ob_zhu']

bigdf = pd.concat([dfdict[d]['df'] for d in dfdict if d not in ignore_datasets])
# Fill NaN's with zeros (i.e. unobserved OTUs)
bigdf = bigdf.fillna(0.0)
bigmeta = pd.concat([dfdict[d]['meta'] for d in dfdict
    if d not in ignore_datasets])

# excludes: 'postFMT_CDI', None, ' '
diseases = ['CIRR', 'MHE', 'HIV', 'T1D', 'EDD', 'CRC', 'ASD', 'CDI',
            'PSA', 'RA', 'OB', 'OB-NASH', 'nonCRC', 'OW',
            'nonCDI', 'PAR', 'CD', 'UC', 'nonNASH-OB', 'NASH']
classes_list = [['H'], diseases]
[h_smpls, dis_smpls] = fio.get_samples(bigmeta, classes_list)

random_state = args.random_state
rf, X, Y = prep_classifier(bigdf, h_smpls, dis_smpls, random_state)

results = cv_and_roc(rf, X, Y, num_cv=args.n_cv, random_state=random_state)

all_results = []

resultsdf = pd.DataFrame(data=np.array((results['mean_fpr'],
                                        results['mean_tpr'])).T,
                         columns=['mean_fpr', 'mean_tpr'])
resultsdf['roc_auc'] = results['roc_auc']
resultsdf['fisher_p'] = results['fisher_p']
resultsdf['dataset'] = 'overall'
resultsdf['H_smpls'] = len(h_smpls)
resultsdf['dis_smpls'] = len(dis_smpls)
resultsdf['num_features'] = bigdf.shape[1]
all_results.append(resultsdf)

## Extract individual predictions
preds = pd.DataFrame(columns=['y_true', 'y_prob'])
preds['y_true'] = results['y_true']
preds['y_prob'] = results['y_prob']
preds['sample'] = h_smpls + dis_smpls
preds = pd.merge(preds, bigmeta[['DiseaseState', 'dataset']],
                 left_on='sample', right_index=True)
preds['disease_class'] = preds['dataset'].map(lambda x: x.split('_')[0])
preds['disease_class'] = preds['disease_class'].replace('edd', 'cdi')

# Calculate disease-wise AUC from the full classifier
for g, subdf in preds.groupby('dataset'):
    if g == 'nash_zhu':
        h = subdf.query('DiseaseState == "H"')
        nash = pd.concat((subdf.query('DiseaseState == "NASH"'), h))
        fpr_nash, tpr_nash, _ = roc_curve(nash['y_true'], nash['y_prob'])

        tmp = pd.DataFrame(data=np.array((fpr_nash, tpr_nash)).T,
                           columns=['fpr', 'tpr'])
        tmp['dataset'] = 'nash_zhu'
        tmp['auc'] = auc(fpr_nash, tpr_nash)
        tmp['H_smpls'] = nash.query('y_true == 0').shape[0]
        tmp['dis_smpls'] = nash.query('y_true == 1').shape[0]
        all_results.append(tmp)

        ob = pd.concat((subdf.query('DiseaseState == "OB"'), h))
        fpr_ob, tpr_ob, _ = roc_curve(ob['y_true'], ob['y_prob'])

        tmp = pd.DataFrame(data=np.array((fpr_ob, tpr_ob)).T,
                           columns=['fpr', 'tpr'])
        tmp['dataset'] = 'ob_zhu'
        tmp['auc'] = auc(fpr_ob, tpr_ob)
        tmp['H_smpls'] = ob.query('y_true == 0').shape[0]
        tmp['dis_smpls'] = ob.query('y_true == 1').shape[0]
        all_results.append(tmp)
    else:
        fpr, tpr, _ = roc_curve(subdf['y_true'], subdf['y_prob'])
        roc_auc = auc(fpr, tpr)
        tmp = pd.DataFrame(data=np.array((fpr, tpr)).T,
                           columns=['fpr', 'tpr'])
        tmp['dataset'] = g
        tmp['auc'] = roc_auc
        tmp['H_smpls'] = subdf.query('y_true == 0').shape[0]
        tmp['dis_smpls'] = subdf.query('y_true == 1').shape[0]
        all_results.append(tmp)

all_results_df = pd.concat(all_results)
all_results_df.to_csv(args.out_file, sep='\t', index=False)
