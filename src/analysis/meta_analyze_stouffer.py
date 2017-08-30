#!/usr/bin/env python
"""
Convert individual q-values into an overall q-value using
Stouffer's method, and return the "core" bacteria.
"""

import os
import argparse
import pandas as pd
import numpy as np
import copy
import scipy.stats as stats

def pvals_to_long(pvals):
    """
    Given dataframe with signed p-values, convert to longform
    with columns: otu, study, direction, pval, sample_size.

    Parameters
    ----------
    pvals : pandas DataFrame
        Genera in rows, studies in columns, signed pvalues in values.
        Positive indicates higher in disease, negatives is higher in healthy.

    Returns
    -------
    longpvals : pandas DataFrame
        Tidy dataframe with columns otu, study, direction, and pval (for that
        direction)
    """
    pvals.index.name = 'otu'
    pvals = pvals.reset_index()
    longpvals = pd.melt(pvals, id_vars='otu', var_name='dataset',
                        value_name='signed_qvalue').dropna()

    # Convert all p-values to health-associated pvalue
    # Original p-values were calculated from KW test, making them two-sided.
    # If the pvalue is negative, then abs(p)/2 is the health-associated pval.
    # If the pvalue is positive, then 1 - abs(p)/2 is the health-associated
    # pvalue.
    p_to_healthy = lambda x: abs(x)/2.0 if x <= 0  else 1-abs(x)/2.0
    longpvals['q'] = longpvals['signed_qvalue'].map(p_to_healthy)
    longpvals['direction'] = 'healthy'

    # Now add the disease-associated qvalues
    disqs = copy.deepcopy(longpvals)
    disqs['direction'] = 'disease'
    disqs['q'] = 1 - disqs['q']

    longpvals = pd.concat((longpvals, disqs))

    return longpvals

parser = argparse.ArgumentParser()
parser.add_argument('qvalues', help='file with qvalues; genera in rows, '
     + 'datasets in columns')
parser.add_argument('dataset_info', help='file with sample sizes; '
    + ' datasets in rows, at least column "total" with sample size.')
parser.add_argument('combined_qvalues', help='out file to write combined '
    + 'qvalues to.')
parser.add_argument('core_bugs', help='out file with +/- 1 indication '
    + 'of core genera.')
parser.add_argument('--pthresh', help='significance threshold [default: '
    + '%(default)s]', default=0.05, type=float)
parser.add_argument('--exclude-nonhealthy', help='flag to exclude '
    + 'studies without healthy controls and hiv_lozupone from the '
    + 'overall cross-disease meta-analysis', action='store_true')
args = parser.parse_args()

pvals = pd.read_csv(args.qvalues, sep='\t', index_col=0)
dataset_info = pd.read_csv(args.dataset_info, sep='\t')

if args.exclude_nonhealthy:
    to_exclude = ['ibd_papa', 'ibd_gevers', 'hiv_lozupone']
    pvals = pvals.drop(to_exclude, axis=1)

## Convert wide pvals df (genera X datasets) into tidy df
longpvals = pvals_to_long(pvals)
# Add sample sizes
longpvals = pd.merge(longpvals, dataset_info[['dataset', 'total']])

## For each [otu, direction] group, combine one-tailed qvalues
## with Stouffer's method
metap = []
for i, subdf in longpvals.groupby(['otu', 'direction']):
    if subdf.shape[0] > 1:
        otu = i[0]
        direction = i[1]
        numstudies = subdf.shape[0]
        studies_str = ','.join(subdf['dataset'])
        # Stouffer's weighted z-score test, weight by sqrt(sample_size)
        z, p = stats.combine_pvalues(subdf['q'].astype(float),
            method='stouffer', weights=subdf['total'].apply(np.sqrt))
        metap.append([otu, direction, z, p, numstudies, studies_str])
metap = pd.DataFrame(metap,
    columns=['otu', 'direction', 'z', 'combined_p', 'num_studies', 'studies'])
metap.to_csv(args.combined_qvalues, sep='\t', index=False)

## Using qthresh, convert combined pvalues to +/- 1 and write to file
def p_to_sig(row, pthresh):
    """
    Convert row with columns 'combined_p' and 'direction' into +/- 1 indicator.
    """
    dirs = {'healthy': -1, 'disease': 1}
    if row['combined_p'] < pthresh:
        return dirs[row['direction']]
    else:
        return np.nan

metap['overall'] = metap.apply(lambda row: p_to_sig(row, args.pthresh), axis=1)

# Get rid of duplicates (bc each OTU has a 'healthy' and 'disease' entry...)
# Make a dataframe with only 'otu' and 'overall' columns
overall = pd.merge(
    pd.DataFrame(metap['otu'].drop_duplicates()),
    pd.DataFrame(metap[['overall', 'otu']].dropna(subset=['overall'])),
    how='outer')
overall.to_csv(args.core_bugs, sep='\t', index=False)
