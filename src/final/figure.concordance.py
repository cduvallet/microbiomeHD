#!/usr/bin/env python
"""
This script plots the pairwise p-values for concordance
"""
import argparse
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
import Formatting as fmt


def get_pvals_df(df):
    """
    Calculate pvalues from the dataframe with the empirical distribution.

    Parameters
    ----------
    df : pandas DataFrame
        has columns ['study1', 'study2', 'qthresh', 'metric', 'n_concord',
        'dis1', 'dis2']. (Produced from concordance_analysis.py)
        The 'metric' column indicates whether the value in 'n_concord' is the
        actual value ("observed") or a shuffled, expected value ("expected").

    Returns
    -------
    pvalsdf : pandas DataFrame
        has columns ['study1', 'study2', 'p'] where 'p' is the p value of the
        expected n_concord being >= the observed n_concord.
        Note that pvalsdf contains *redundant* study1-study2 pairs (so that
        calling pvalsdf.pivot() yields a symmetric matrix)
    """

    # Make pvalues dataframe from the empirical distribution
    pvals = []
    for g, subdf in df.groupby(['study1','study2', 'qthresh']):
        n_obs = subdf.query('metric == "observed"')['n_concord'].values
        n_exp = subdf.query('metric == "expected"')['n_concord'].values
        p = sum(n_exp >= n_obs)/float(len(n_exp))

        pvals.append(
            list(g) + [p, subdf['dis1'].iloc[0], subdf['dis2'].iloc[0]])
    pvalsdf = pd.DataFrame(data=pvals,
        columns=['study1', 'study2', 'qthresh', 'p', 'dis1', 'dis2'])
    # Add study1-study2 flipped pairs so pivoted df is symmetric
    pvalsdf = pd.concat(
        (pvalsdf,
         pvalsdf[['study1', 'study2', 'p']].rename(
            columns={'study1': 'study2', 'study2': 'study1'})))\
        .drop_duplicates(['study1', 'study2', 'p'])
        # the within-disease comparisons get doubled here
    return pvalsdf

## Arguments
p = argparse.ArgumentParser()
p.add_argument('concordance', help='tidy file with concordance results')
p.add_argument('dataset_info', help='file with dataset info, for use in '
    + 'ordering datasets for plotting (in Formatting.py)')
p.add_argument('fout', help='path to write figure to')
args = p.parse_args()

## Read in stuff
df = pd.read_csv(args.concordance, sep='\t')
dataset_info = pd.read_csv(args.dataset_info, sep='\t')
pvalsdf = get_pvals_df(df)

## Prepare for plotting
# Get dataset orders, as in Figure 1
dataset_info = dataset_info\
    .replace('edd_singh', 'cdi_singh')\
    .replace('noncdi_schubert', 'cdi_schubert2')
_, dataset_order = fmt.get_dataset_order(dataset_info)
# Set up matrix to plot
toplot = pvalsdf\
    .replace('edd_singh', 'cdi_singh')\
    .replace('noncdi_schubert', 'cdi_schubert2')\
    .pivot('study2', 'study1', 'p')
toplot = toplot.loc[dataset_order, dataset_order]
toplot = pd.DataFrame(data=np.tril(np.log10(toplot.values + 1e-4)),
    columns=toplot.columns, index=toplot.index)
toplot = toplot.replace(0, np.nan)
# Get labels for each dataset
labeldict = fmt.get_labeldict(toplot.index)
vmin = -4 #should match the pseudo-zero addition
vmax = 0

## Set up grid spec - will make boxes around each disease
sns.set_style('white')
fig = plt.figure(figsize=(7,6))
gs = gridspec.GridSpec(toplot.shape[0], toplot.shape[1],
    left=0.25, right=0.95, bottom=0.2, wspace=0.2, hspace=0.2)

# Get all the diseases
alldis = []
for d in dataset_order:
    d = d.split('_')[0]
    if d not in alldis:
        alldis.append(d)

for dis1 in alldis:
    for dis2 in alldis[alldis.index(dis1):]:

        yindices = [i for i in range(len(dataset_order))
                    if dataset_order[i].startswith(dis1)]
        xindices = [i for i in range(len(dataset_order))
                    if dataset_order[i].startswith(dis2)]

        # Size each axis (i.e. number of gridspec slots) according to the
        # number of studies in that disease
        ax = plt.subplot(gs[xindices[0]:xindices[0]+len(xindices),
                            yindices[0]:yindices[0]+len(yindices)])

        cols = [i for i in toplot.index if i.startswith(dis1)]
        rows = [i for i in toplot.columns if i.startswith(dis2)]
        subdf = toplot.loc[rows, cols]
        _ = ax.imshow(subdf.values, interpolation='nearest', aspect='auto',
                      cmap=sns.light_palette("gray", reverse=True,
                                             as_cmap=True),
                      vmin=vmin, vmax=vmax)

        #ax.spines['bottom'].set_color('gray')
        #ax.spines['left'].set_color('gray')
        #ax.spines['right'].set_color('gray')
        #ax.spines['top'].set_color('gray')
        sns.despine(ax=ax, left=True, bottom=True)

        # Labels
        if dis1 == 'cdi':
            if dis2 == 'cdi':
                labels = [''] + [labeldict[i] for i in subdf.index[1:]]
            else:
                labels = [labeldict[i] for i in subdf.index]
            ax.set_yticks(range(0, subdf.shape[0]))
            _ = ax.set_yticklabels(labels, fontsize=8, va='center')
        else:
            ax.set_yticklabels([])

        if dis2 == 'par' and dis1 != 'par':
            labels = [labeldict[i] for i in subdf.columns]
            ax.set_xticks(range(0, subdf.shape[1]))
            # Align rotated labels with actual bar
            ax.set_xticks(ax.get_xticks() + 0.4)
            _ = ax.set_xticklabels(labels, rotation=45, fontsize=8, ha='right')
        else:
            ax.set_xticklabels([])

# Save
fig.savefig(args.fout)

## Colorbar
fig, ax = plt.subplots(figsize=(5,1))
_ = ax.imshow(np.tile(np.linspace(-4, 0, num=100), (5, 1)), interpolation='nearest', aspect='auto', cmap=sns.light_palette("gray", reverse=True, as_cmap=True))
sns.despine(ax=ax, left=True, bottom=True)
#ax.set_xticks([0, 100])
#ax.set_xticklabels([r'$p < 0.0001', r'$p = 1$'], rotation=90)
ax.set_xticklabels([])
ax.set_yticklabels([])
fig.savefig('final/figures/concordance_colorbar.pdf')
