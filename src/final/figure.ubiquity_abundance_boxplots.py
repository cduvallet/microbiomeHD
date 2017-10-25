#!/usr/bin/env python
"""
This script plots boxplots of abundance/ubiquity for each of the genera,
separated by core status.
"""
import argparse

import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('TKAgg')
import seaborn as sns
import matplotlib.pyplot as plt


def plot_ubiq_abun_boxplot(tidy, metric, calculation):
    """
    Plot boxplot where x-axis is 'overall_significance' of genus, and values
    are either ubiquity or abundance in tidy (with the respective metric and
    calculation type)

    Parameters
    ----------
    tidy : pandas dataframe
        has columns overall_significance, value, patient, metric, and calculation
    metric : str
        'abundance' or 'ubiquity'
    calculation: str
        'from_pooled_mean' or 'mean_of_datasets'

    Returns
    -------
    ax : Axis object
    """
    fig, ax = plt.subplots(figsize=(5.5,4))
    tmp = tidy.query('metric == @metric')\
              .query('calculation == @calculation')\
              .query('patient == "total"')

    boxprops = {'edgecolor': 'k', 'facecolor': 'w'}
    lineprops = {'color': 'k'}

    # Plot log10(abundance)
    if metric == 'abundance':
        tmp.loc[tmp.index, 'value'] = tmp['value'].apply(np.log10)

    sns.boxplot(data=tmp, x='overall_significance', y='value',
                fliersize=0, ax=ax, color='w',
                order=['health', 'disease', 'mixed', 'not_sig'],
                **{'boxprops': boxprops, 'medianprops': lineprops,
                   'whiskerprops': lineprops, 'capprops': lineprops})

    sns.stripplot(data=tmp, x='overall_significance', y='value',
                  jitter=True, linewidth=0.6, split=True, ax=ax,
                  order=['health', 'disease', 'mixed', 'not_sig'],
                  color='w')
    return fig, ax


p = argparse.ArgumentParser("Plots boxplots of abundance or ubiquity")
p.add_argument('ubiquity', help='path to tidy dataframe with ubiquity values')
p.add_argument('metric', help='which metric to plot',
    choices=['abundance', 'ubiquity'])
p.add_argument('out', help='file to write figure to')
p.add_argument('--calc', help='method used to calculate the metric '
    + '[default: %(default)s]',
    choices=['from_pooled_mean', 'mean_of_datasets'],
    default='mean_of_datasets')
args = p.parse_args()

tidy = pd.read_csv(args.ubiquity, sep='\t')

sns.set_style('white')
fig, ax = plot_ubiq_abun_boxplot(tidy, args.metric, args.calc)
ax.set_title('')
ax.set_xticklabels(['Non-specific\nhealth', 'Non-specific\ndisease',
                    'Non-specific\nmixed', 'Not\nnon-specific'],
                    fontsize='medium')
ax.set_xlabel('')#, fontsize='x-large')

if args.metric == 'abundance':
    ax.set_ylabel('log10(Abundance)', fontsize='x-large')
else:
    ax.set_ylabel('Ubiquity', fontsize='x-large')

if args.metric == 'ubiquity':
    ax.set_ylim([-0.05, 1.05])

ax.set_yticklabels(ax.get_yticks(), fontsize='large')
fig.tight_layout()
fig.savefig(args.out)
