#!/usr/bin/env python
"""
This script runs genus-level classifiers on all datasets.
It writes a tidy data file with the points to plot each dataset's
ROC curve and also the corresponding AUC for each.
"""

import argparse
import pandas as pd
import numpy as np
from sklearn.metrics import cohen_kappa_score

import os, sys
# Add src/util to path and import modules from files there
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
import FileIO as fio
from util import collapse_taxonomic_contents_df, prep_classifier, cv_and_roc

def results2df(results, dataset, n_ctrl, n_case, n_features):
    """
    Converts the results dictionary into tidy dataframe format.

    Parameters
    ----------
    results : dict
        results from cv_and_roc function
    dataset : str
        dataset ID
    n_ctrl, n_case, n_features : int
        number of controls, cases, and features

    Returns
    -------
    resultsdf : pandas DataFrame
        dataframe with 'mean_fpr', 'mean_tpr', 'fisher_p', and 'roc_auc' columns
        from the results dict, 'kappa' from
        cohen_kappa_score(results['y_preds']), and 'dataset', 'H_smpls',
        'dis_smpls', and 'num_features' from the input parameters
    """
    # Directly calling pd.DataFrame.from_dict doesn't work because
    # this dictionary contains arrays, matrices, etc..
    resultsdf = pd.DataFrame(data=np.array((results['mean_fpr'],
                                            results['mean_tpr'])).T,
                             columns=['mean_fpr', 'mean_tpr'])
    resultsdf['roc_auc'] = results['roc_auc']
    resultsdf['fisher_p'] = results['fisher_p']

    resultsdf['dataset'] = dataset
    resultsdf['H_smpls'] = n_ctrl
    resultsdf['dis_smpls'] = n_case
    resultsdf['num_features'] = n_features

    resultsdf['kappa'] = cohen_kappa_score(
        results['y_true'], results['y_preds'])

    return resultsdf

p = argparse.ArgumentParser()
p.add_argument('datadir', help='directory with clean OTU tables and metadata.')
p.add_argument('outfile', help='out file with tidy RF results')
p.add_argument('--randomstate', help='random state seed. [default: '
    + '%(default)s]', default=12345)
p.add_argument('--core', help='path to file with core bugs, which has genera '
    + 'in rows and column "overall" with non-NaN values indicating coreness. '
    + 'If given, builds classifiers with only core bugs as features.',
    default=None)
p.add_argument('--subset', help='file with list of dataset IDs that '
    + 'should be analyzed, if you want to analyze only a subset of the '
    + 'datasets in clean_data_dir. Dataset IDs should match IDs in file '
    + 'names, and be one per line.', default=None)
p.add_argument('--split-cases', help='flag to analyze each case type '
    + 'separately.', action='store_true')
args = p.parse_args()

dfdict = fio.read_dfdict_data(args.datadir, subset=args.subset)

## If we should make classifiers using only the core bugs, grab those now
if args.core is not None:
    overall = pd.read_csv(args.core, sep='\t', index_col=0)
    core_bugs = overall.dropna().index.tolist()

tidyresults = []

## Classify each dataset
print('Classifying datasets...')
for dataset in dfdict.keys():
    df = dfdict[dataset]['df']
    meta = dfdict[dataset]['meta']

    # Prepare OTU table: collapse to genus, keep only core (if required)
    df = collapse_taxonomic_contents_df(df,'genus')
    if args.core is not None:
        keep_bugs = [i for i in core_bugs if i in df.columns]
        df = df[keep_bugs]

    if args.split_cases:
        # Get samples
        classes_list = fio.get_classes(meta, dataset)
        # Iterate through case groups
        for dis_label in classes_list[1]:
            # e.g. old dataset = ibd_alm, new dataset = uc_alm
            newdataset = dis_label.lower() + '_' + dataset.split('_')[1]
            print(newdataset),

            # Get samples
            sub_list = [classes_list[0], [dis_label]]
            H_smpls, dis_smpls = fio.get_samples(meta, sub_list)

            # Make RF
            rf, X, Y = prep_classifier(df, H_smpls, dis_smpls, args.randomstate)
            results = cv_and_roc(rf, X, Y, random_state=args.randomstate)

            # Update results
            resultsdf = results2df(results, newdataset,
                                   len(H_smpls), len(dis_smpls), df.shape[1])
            tidyresults.append(resultsdf)

    else:
        print(dataset),
        classes_list = fio.get_classes(meta, dataset)
        H_smpls, dis_smpls = fio.get_samples(meta, classes_list)

        rf, X, Y = prep_classifier(df, H_smpls, dis_smpls, args.randomstate)
        results = cv_and_roc(rf, X, Y, random_state=args.randomstate)

        resultsdf = results2df(results, dataset,
                               len(H_smpls), len(dis_smpls), df.shape[1])

        tidyresults.append(resultsdf)

tidydf = pd.concat(tidyresults, ignore_index=True)
tidydf.to_csv(args.outfile, sep='\t', index=False)
