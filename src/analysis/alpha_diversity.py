#!/usr/bin/env python
"""
This script calculates the alpha diversity for each non-collapsed
OTU table.
"""
import numpy as np
import pandas as pd
import argparse

import skbio.diversity.alpha as alph

from scipy.stats import ranksums, ttest_ind
from scipy.stats.mstats import kruskalwallis
from statsmodels.sandbox.stats.multicomp import multipletests

import os, sys
# Add src/util to path and import modules from files there
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
from FileIO import get_dataset_ids, read_dataset_files
from util import raw2abun


def get_all_pvals(df, groupcol, valuecol, method='kruskalwallis'):
    """
    Returns pairwise p-values between all groups in the column `groupcol`.

    Parameters
    ----------
    df : pandas dataframe
        tidy dataframe with labels in `groupcol` and values in `valuecol`
    groupcol, valuecol : str
        columns in df
    method : str {'kruskalwallis', 'ranksums', 'wilcoxon', 'ttest_ind'}
        statistical method for comparison. Default is 'kruskalwallis'

    Returns
    -------
    pvals : dict
        dictionary with 'group1_vs_group2' as the keys and p-value as the values
    """

    pvals = {}

    ## Get all pairwise combinations
    grps = list(set(df[groupcol]))
    for g1 in grps:
        for g2 in grps[grps.index(g1)+1:]:
            if g1 != g2:
                ## Grab values
                x = df[df[groupcol] == g1][valuecol]
                y = df[df[groupcol] == g2][valuecol]

                ## Calculate p value
                if method == 'ranksums' or method == 'wilcoxon':
                    pfun = ranksums
                elif method == 'ttest_ind':
                    pfun = ttest_ind
                else:
                    pfun = kruskalwallis
                try:
                    _, p = pfun(x, y)
                except:
                    # Should probably have better error handling here...
                    p = np.nan

                ## Store p value
                pvals[g1 + '_vs_' + g2] = p
    return pvals

def get_layered_pvals(df, groupcol, valuecol, subset_by,
                      pval_method='kruskalwallis'):
    """
    Get pvalues for all pairwise combinations in groupcol.
    Performs calculating separately for each group in subset_by columns.
    In other words, this is a wrapper for groupby(subset_by) + get_all_pvals().

    Parameters
    ----------
    df : pandas dataframe
        tidy dataframe with labels in `groupcol` and values in `valuecol`
    groupcol, valuecol : str
        columns in df
    subset_by : str
        column to group by
    pval_method : str {'kruskalwallis', 'ranksums', 'wilcoxon', 'ttest_ind'}
        statistical method for comparison. Default is 'kruskalwallis'

    Returns
    -------
    pvals : dict
        multi-level dictionary, with outside keys as the unique values in
        df[subset_by] and the inner values as in get_all_pvals()
    """

    pvals = {}
    for s, subdf in df.groupby(subset_by):
        pvals[s] = get_all_pvals(subdf, groupcol, valuecol,
                                      method=pval_method)
    return pvals

def alpha_diversity(df, metric='shannon'):
    """
    Calculate shannon diversity of all samples in df.

    Parameters
    ----------
    df : pandas DataFarme
        dataframe with samples in rows, OTUs in columns
    metric : str
        'shannon', 'chao1', 'simpson'

    Returns
    -------
    alpha : pandas Series
        pandas series with samples in rows
    """

    alphafun = alph.shannon
    if metric == 'chao1':
        alphafun = alph.chao1
    elif metric == 'simpson':
        alphafun = alph.simpson
    elif metric != 'shannon':
        print('Unknown alpha diversity metric. Doing Shannon Index.')

    return df.apply(alphafun, axis=1)

def make_alpha_df(df, meta, study, metric):
    """
    Make the alpha diversity tidy dataframe.

    Parameters
    ----------
    df: pandas DataFrame
        samples in rows, OTUs in columns
    meta : pandas DataFrame
        samples in rows, 'DiseaseState' column
    study : str
        study label
    metric : str
        alpha diversity metric

    Returns
    -------
    alpha : pandas DataFrame
        tidy dataframe with ['sample', 'alpha', 'alpha_metric', 'study',
        'DiseaseState'] columns
    """
    alpha = alpha_diversity(df, metric).reset_index()
    alpha.columns = ['sample', 'alpha']
    alpha['alpha_metric'] = metric
    alpha['study'] = dataset
    # Label diseases based on index, which is now in first column
    alpha['DiseaseState'] = alpha['sample']\
        .apply(lambda x: meta.loc[x, 'DiseaseState'])
    return alpha


p = argparse.ArgumentParser()
p.add_argument('datadir', help='directory with clean OTU tables and metadata.')
p.add_argument('alphas_out', help='out file with all alpha diversities for '
    + 'all samples')
p.add_argument('pvals_out', help='out file with all alpha diversity p-values')

args = p.parse_args()

## Data in datadir has already been cleaned up, but not collapsed to genus level
datadir = args.datadir

datasetids = get_dataset_ids(datadir)
alphas = []

for dataset in datasetids:
    print(dataset),
    ## Read dataset
    df, meta = read_dataset_files(dataset, datadir)

    for metric in ['shannon', 'chao1', 'simpson']:
        alpha = make_alpha_df(df, meta, dataset, metric)
        alphas.append(alpha)

alphasdf = pd.concat(alphas, ignore_index=True)

alphasdf.to_csv(args.alphas_out, sep='\t', index=False)

# Because I'm using the entire OTU table, some of these samples don't have disease metadata.
# I don't want to compare "NaN" labeled samples with anything because they mean nothing.
alphasdf = alphasdf.query('DiseaseState != " "').dropna(subset=['DiseaseState'])

pvals = []
for g, subdf in alphasdf.groupby('alpha_metric'):
    pval = get_layered_pvals(subdf, 'DiseaseState', 'alpha', 'study')
    pval = pd.DataFrame.from_dict(pval).stack().reset_index()
    pval.columns = ['comparison', 'study', 'p']
    pval['q'] = multipletests(pval['p'])[1]
    pval['alpha_metric'] = g
    pvals.append(pval)

pvalsdf = pd.concat(pvals)
pvalsdf.to_csv(args.pvals_out, sep='\t', index=False)
