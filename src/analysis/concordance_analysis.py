#!/usr/bin/env python
"""
This script looks at the expected concordance between all studies in
all pairwise disease comparisons.
"""
import argparse

import pandas as pd
import numpy as np

from scipy.stats import fisher_exact, spearmanr, kendalltau
from sklearn.metrics import cohen_kappa_score

# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
from util import shuffle_col

def empirical_pval(series1, series2, nreps=1000):
    """
    Get the empirical pvalue of n_concordance between values in series1
    and series2.

    Parameters
    ----------
    series1, series2 : pandas Series
        NaN structure is respected during shuffling.
    nreps : int
        number of reps to build null from

    Returns
    -------
    effect : float
        n_obs - mean(n_expected)
    p : float
        sum(n_expected >= n_observed) where n_expected is the null distribution
        of the number of concordant values
    """

    observed = sum(series1 == series2)

    dist = []
    for i in range(nreps):
        # shuffle both series
        shuffled1 = shuffle_col(series1)
        shuffled2 = shuffle_col(series2)
        # Merge the shuffled series so they have the same indices
        shuffled = pd.concat((shuffled1, shuffled2), axis=1)
        dist.append(sum(shuffled.iloc[:, 0] == shuffled.iloc[:, 1]))

    effect = observed - np.mean(dist)
    p = sum(dist >= observed)/float(len(dist))

    return effect, p

def concordance(series1, series2, method, nreps=1000):
    """
    Measures the concordance between two pandas Series and returns a pvalue
    and measure of concordance.

    Parameters
    ----------
    series1, series2 : pandas Series
        Series with matching indexes.
    method : str
        ['fisher', 'spearman', 'kendalltau', 'empirical', 'cohen']
    nreps : int
        number of repititions to build the null. Only needed if method is
        'empirical'

    Returns
    -------
    measure : float
        some sort of measure of concordance (e.g. r for the correlation
        methods, n_observed - mean(n_expected) for empirical, etc)
    p : float
        p value of observed concordance between series1 and series2
    """

    if method == 'fisher':
        # Note: this automatically ignores any bugs which were not present
        # in both series.
        mat = pd.crosstab(series1, series2)
        return fisher_exact(mat)

    elif method == 'spearman':
        return spearmanr(series1, series2)

    elif method == 'kendalltau':
        return kendalltau(series1, series2, nan_policy='omit')

    elif method == 'empirical':
        return empirical_pval(series1, series2, nreps)

    elif method == 'cohen':
        tmp = pd.concat((series1, series2), axis=1).dropna()
        return cohen_kappa_score(tmp.iloc[:, 0], tmp.iloc[:, 1]), np.nan

    else:
        raise ValueError('Unknown concordance method.')

p = argparse.ArgumentParser()
p.add_argument('qvals', help='path to file with signed qvalues. Datasets '
    + 'in columns, genera in rows, signed qvalues in values.')
p.add_argument('--nreps', help='number of shuffles to build empirical null. '
    + '[default: %(default)s]',
    default=1000, type=int)
p.add_argument('fout', help='file to write pvalues to.')
args = p.parse_args()

# Read in qvalues
df = pd.read_csv(args.qvals, sep='\t', index_col=0)
# Rename edd_singh to cdi_singh
df = df.rename(columns={'edd_singh': 'cdi_singh',
                        'noncdi_schubert': 'cdi_schubert2'})

# Convert df to effect directions
# TODO: include different thresholds? Perhaps as an argument to fxn?
df = np.sign(df).replace(0, np.nan)

methods = ['fisher', 'spearman', 'kendalltau', 'empirical', 'cohen']

# Iterate over all pairwise studies
allstudies = list(df.columns)
results = []
for study1 in allstudies:
    print(study1)
    for study2 in allstudies[allstudies.index(study1)+1:]:
        dis1 = study1.split('_')[0]
        dis2 = study2.split('_')[0]

        series1 = df[study1]
        series2 = df[study2]

        for method in methods:
            measure, p = concordance(series1, series2, method, args.nreps)
            results.append([dis1, dis2, study1, study2, measure, p, method])

resultsdf = pd.DataFrame(data=results,
    columns=['dis1', 'dis2', 'study1', 'study2', 'measure', 'p', 'method'])
resultsdf.to_csv(args.fout, sep='\t', index=False)
