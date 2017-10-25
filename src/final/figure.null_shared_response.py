#!/usr/bin/env python
"""
This script plots the empirical null number of non-specific responders
for each n_diseases heuristic.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import argparse

def plot_null(fnull, fcore, ax=None):
    """
    Plot the null shared response overlayed with the actual number of
    shared responders.
    """
    null = pd.read_csv(fnull, sep='\t')
    core = pd.read_csv(fcore, sep='\t', index_col=0)
    ncore = core.groupby('overall').size()

    for i in [-1, 0, 1]:
        if i not in ncore:
            ncore.loc[i] = 0

    order = ['health', 'mixed', 'disease']

    if ax is None:
        fig, ax = plt.subplots()

    sns.stripplot(data=null, x='type', y='n', order=order, jitter=True, ax=ax, alpha=0.2)

    ax.scatter([0, 1, 2], [ncore.loc[-1], ncore.loc[0], ncore.loc[1]],
               c='k', marker='D', s=25, zorder=10)

    return ax

p = argparse.ArgumentParser()
p.add_argument('null_stem', help='filepath stem for file with "rep", "type", '
    + 'and "n" columns, where "n" is the number of observed non-specific '
    + 'responders after shuffling qvalues. Full files will have format '
    + ' null_stem.n_diseases.txt where n is 2, 3, 4, or 5.')
p.add_argument('core_stem', help='filepath stem for file with "overall" '
    + 'column, containing +/- 1 or 0 values indicating nonspecific status. '
    + 'Full files will have format '
    + 'core_stem.2_diseases.across_all_diseases.txt, where n is 2, 3, 4, or 5.')
p.add_argument('out', help='out file to save figure to')

args = p.parse_args()

# To help with printing pvalues later
d = {'health': -1, 'mixed': 0, 'disease': 1}
order = ['health', 'mixed', 'disease']

sns.set_style('white')
fig, ax = plt.subplots(2,2, figsize=(5,5))
ax = ax.flatten()
n_diseases = [2,3,4,5]
for i in range(4):
    n = n_diseases[i]
    fnull = '{}.{}_diseases.txt'.format(args.null_stem, n)
    fcore ='{}.{}_diseases.across_all_diseases.txt'.format(args.core_stem, n)
    ax[i] = plot_null(fnull, fcore, ax=ax[i])
    ax[i].set_title('{} diseases'.format(n))
    # x label
    if i in [0, 1]:
        ax[i].set_xticklabels([])
    ax[i].set_xlabel('')
    # y label
    if i in [1, 3]:
        ax[i].set_ylabel('')
    else:
        ax[i].set_ylabel('N non-specific')
    # y limits
    newymin = ax[i].get_yticks()[0]/5.0
    oldymax = ax[i].get_yticks()[-1]
    ytickrange = ax[i].get_yticks()[-1] - ax[i].get_yticks()[-2]
    newymax = oldymax + ytickrange
    ax[i].set_ylim([newymin, newymax])
    if i in [2]:
        ax[i].set_yticks([0, 1, 2, 3, 4])
    if i in [3]:
        ax[i].set_yticks([0, 1])

    # print pvalues
    null = pd.read_csv(fnull, sep='\t')
    core = pd.read_csv(fcore, sep='\t', index_col=0)
    ncore = core.groupby('overall').size()
    for j in [-1, 0, 1]:
        if j not in ncore:
            ncore.loc[j] = 0
    for g in order:
        p = float(sum(null.groupby('type').get_group(g)['n'] >= ncore[d[g]])) / len(null.groupby('type').get_group(g)['n'])
        ax[i].text(s='{:.2f}'.format(p), x=(d[g]+1),
            y=(oldymax + 0.25*ytickrange), ha='center', fontsize=8)

sns.despine(fig=fig)
fig.tight_layout()
fig.savefig(args.out)
