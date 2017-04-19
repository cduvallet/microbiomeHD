#!/usr/bin/env python
"""
This script runs genus-level classifiers on all datasets.
It writes a tidy data file with the points to plot each dataset's
ROC curve and also the corresponding AUC for each.
"""

import argparse
import pandas as pd
import numpy as np

import os, sys
# Add src/util to path and import modules from files there
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
from FileIO import read_dfdict_data
from util import collapse_taxonomic_contents_df, prep_classifier, cv_and_roc


p = argparse.ArgumentParser()
p.add_argument('datadir', help='directory with clean OTU tables and metadata.')
p.add_argument('outfile', help='out file with tidy RF results')
p.add_argument('--randomstate', help='random state seed', default=12345)
args = p.parse_args()

dfdict = read_dfdict_data(args.datadir)

tidyresults = []

## Classify each dataset
print('Classifying datasets...')
for dataset in dfdict.keys():
    print(dataset),
    df = dfdict[dataset]['df']
    meta = dfdict[dataset]['meta']
    H_smpls = dfdict[dataset]['H_smpls']
    dis_smpls = dfdict[dataset]['dis_smpls']

    df = collapse_taxonomic_contents_df(df,'genus')

    rf, X, Y = prep_classifier(df, H_smpls, dis_smpls, args.randomstate)
    results = cv_and_roc(rf, X, Y, random_state=args.randomstate)
    # Directly calling pd.DataFrame.from_dict doesn't work because
    # this dictionary contains arrays, matrices, etc..
    resultsdf = pd.DataFrame(data=np.array((results['mean_fpr'],
                                            results['mean_tpr'])).T,
                             columns=['mean_fpr', 'mean_tpr'])
    resultsdf['roc_auc'] = results['roc_auc']
    resultsdf['dataset'] = dataset
    resultsdf['H_smpls'] = len(H_smpls)
    resultsdf['dis_smpls'] = len(dis_smpls)
    tidyresults.append(resultsdf)

tidydf = pd.concat(tidyresults, ignore_index=True)
tidydf.to_csv(args.outfile, sep='\t', index=False)
