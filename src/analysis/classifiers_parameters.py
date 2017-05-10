#!/usr/bin/env python
"""
This script builds classifiers across various RandomForestClassifier
parameters.
"""
import argparse
import multiprocessing
import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score

import os, sys
# Add src/util to path and import modules from files there
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
from FileIO import read_dfdict_data
from util import collapse_taxonomic_contents_df, cv_and_roc, prep_classifier


def run_one_rf((dataset, X, Y, n_est, crit, min_split, min_leaf, random_state)):
    """
    Train and test one classifier.

    Note: parameters are actually passed as one tuple in order for this
    function to be used with multiple processes.
    Also note that I could initialize the RF outside of this function, but
    I want to also return the values used for each parameter so I might as
    well do it in here...?

    Parameters
    ----------
    dataset : str
    X : array-like
        training samples
    Y : list
        true labels
    rf : RandomForestClassifier
        initialized classifier instance with the given parameters
    n_est : int
        RandomForestClassifier n_estimators
    crit : str
        RandomForestClassifier criterion
    min_split : int
        RandomForestClassifier min_samples_split
    min_leaf : int
        RandomForestClassifier min_samples_leaf
    random_state : int
        RandomForestClassifier random_state

    Returns
    -------
    results : list
        [
        dataset (str),
        n_est, crit, min_split, min_leaf,
        roc_auc (float, ROC AUC on interpolated ROC curve),
        fisher_p (float, fisher's exact p-value on confusion matrix),
        auc_precision_recall (float, average precision score),
        score (float, RF oob_score)
        ]

    """
    print(dataset, n_est, crit, min_split, min_leaf),
    rf = RandomForestClassifier(n_estimators=n_est,
                                criterion=crit,
                                min_samples_split=min_split,
                                min_samples_leaf=min_leaf,
                                random_state=random_state)


    # Cross-validated results
    try:
        results = cv_and_roc(rf, X, Y, random_state=random_state)
        # Get precision/recall from the cross-validated resutls
        auc_precision_recall = \
            average_precision_score(results['y_true'], results['y_prob'])
    except:
        results = {'roc_auc': np.nan, 'fisher_p': np.nan}
        auc_precision_recall = np.nan

    # Get oob_score for non-cross validated result
    rf = RandomForestClassifier(n_estimators=n_est, criterion=crit,
                                min_samples_split=min_split,
                                min_samples_leaf=min_leaf,
                                random_state=random_state,
                                oob_score=True)
    try:
        score = rf.fit(X, Y).oob_score_
    except:
        score = np.nan

    print('Finished. (AUC = {:.2f})'.format(results['roc_auc']))

    return [dataset, n_est, crit, min_split, min_leaf,
            results['roc_auc'], results['fisher_p'],
            auc_precision_recall, score]

if __name__ == "__main__":
    # I need to wrap this code in if name == main to use multiprocessing.map
    # and avoid recursion stuff.
    p = argparse.ArgumentParser()
    p.add_argument('datadir', help='directory with clean tables')
    p.add_argument('results_out', help='path to file to write results')
    p.add_argument('--random_state', help='random state seed (default: %(default)s)',
                   default=12345)
    args = p.parse_args()

    random_state = args.random_state
    ## Test out a few different RF parameters to make sure results are
    ## relatively constant
    n_estimators = [1000, 10000]

    # The function to measure the quality of a split. Supported criteria are
    # "gini" for the Gini impurity and "entropy" for the information gain.
    criterion = ['gini', 'entropy']

    # The minimum number of samples required to split an internal node:
    # If int, then consider min_samples_split as the minimum number.
    # If float, then min_samples_split is a percentage and ceil(min_samples_split * n_samples) are the minimum number of samples for each split.
    min_samples_split = [2, 3, 4, 0.1]

    # The minimum number of samples required to be at a leaf node:
    # If int, then consider min_samples_leaf as the minimum number.
    # If float, then min_samples_leaf is a percentage and ceil(min_samples_leaf * n_samples) are the minimum number of samples for each node.
    min_samples_leaf = [1, 2, 3]

    dfdict = read_dfdict_data(args.datadir)
    # collapse all to genus level
    for dataset in dfdict:
        dfdict[dataset]['df'] = \
            collapse_taxonomic_contents_df(dfdict[dataset]['df'], 'genus')

    # Set up list of inputs to run_one_rf
    all_dataset = []
    all_X = []
    all_Y = []
    all_nest = []
    all_crit = []
    all_minsplit = []
    all_minleaf = []
    all_random = []
    for dataset in dfdict.keys():

        print(dataset)
        df = dfdict[dataset]['df']
        meta = dfdict[dataset]['meta']
        H_smpls = dfdict[dataset]['H_smpls']
        dis_smpls = dfdict[dataset]['dis_smpls']

        _, X, Y = prep_classifier(df, H_smpls, dis_smpls, random_state)

        for crit in criterion:
            for min_split in min_samples_split:
                for min_leaf in min_samples_leaf:
                    for n_est in n_estimators:
                        all_dataset.append(dataset)
                        all_X.append(X)
                        all_Y.append(Y)
                        all_nest.append(n_est)
                        all_crit.append(crit)
                        all_minsplit.append(min_split)
                        all_minleaf.append(min_leaf)
                        all_random.append(random_state)

    p = multiprocessing.Pool()
    rf_results = p.map(run_one_rf, zip(all_dataset, all_X, all_Y,
                                       all_nest, all_crit, all_minsplit,
                                       all_minleaf, all_random))

    rf_results_df = pd.DataFrame(rf_results,
                                 columns=['dataset', 'n_estimators',
                                          'criterion', 'min_samples_split',
                                          'min_samples_leaf', 'roc_auc',
                                          'fisher_p', 'auc_prec_recall',
                                          'oob_score'])

    rf_results_df.to_csv(args.results_out, sep='\t', index=False)
