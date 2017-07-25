#!/usr/bin/env python
"""
This script plots the alpha diversity by study.
"""
import argparse

import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import numpy as np
from math import floor, ceil

# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
import Formatting as fmt
from FileIO import get_classes

def plot_alphadf(alphasdf, col_order, labeldict, metric='alpha'):
    """
    Plot faceted alpha diversity.

    Parameters
    ----------
    alphasdf : pandas DataFrame
        columns ['study', 'alpha', 'DiseaseState']
    col_order : list
        dataset IDs in the order they should be plotted
    labeldict : dict
        dictionary with {dataset: label}
    mteric : str
        alpha diversity metric, to use in labeling y axis

    Returns
    -------
    fig : Figure
    """
    sns.set_style('white')
    g = sns.FacetGrid(alphasdf, col='study', col_wrap=6,
                      col_order=col_order, sharex=False, sharey=False)
    g = g.map(sns.boxplot, "DiseaseState", "alpha")
    g = g.map(sns.stripplot, "DiseaseState", "alpha", split=True, jitter=True,
              size=5, linewidth=0.6)

    fig = plt.gcf()
    fig.set_size_inches(14.2, 9)

    # Fix y-axis gridlines
    axs = g.axes
    for i in range(len(axs)):
        ax = axs[i]
        yticks = ax.get_yticks()
        # If bottom limit is between 0 and 1 (i.e. not simpson)
        if not (yticks[0] < 1 and yticks[0] > 0):
            ax.set_ylim(floor(yticks[0]), floor(yticks[-1]))
        if yticks[0] < 0:
            ax.set_ylim(0, floor(yticks[-1]))

        yticks = ax.get_yticks()
        if (yticks[0] < 1 and yticks[0] > 0):
            ax.set_yticks(yticks[1::2])
        else:
            ax.set_yticks(yticks[::2])
            # Need some space on the y-axis for p-values
            ax.set_ylim(ax.get_ylim()[0], 1.2*ax.get_ylim()[1])
        # Update title
        oldtitle = ax.get_title()
        newtitle = labeldict[oldtitle.split('=')[1].strip()]
        ax.set_title(newtitle)

        # Update y label
        if i % 6 == 0:
            ax.set_ylabel(metric)

    plt.tight_layout()
    return fig

p = argparse.ArgumentParser()
p.add_argument('alpha', help='file with tidy alpha diversity data')
p.add_argument('out', help='file prefix for alpha diversity figure')
args = p.parse_args()

## Prepare alpha diversity dataframe
alphasdf = pd.read_csv(args.alpha, sep='\t')

# Remove un-annotated samples
alphasdf = alphasdf.replace(' ', np.nan).replace('IBDundef', np.nan).dropna()

# Get rid of noncdi_schubert study and replace cdi_schubert's ignore-nonCDI
alphasdf = alphasdf.query('study != "noncdi_schubert"')
alphasdf = alphasdf.replace('ignore-nonCDI', 'nonCDI')

# Sort datasets alphabetically, with edd_singh grouped with the CDI studies
alphasdf['study'] = alphasdf['study'].replace('edd_singh', 'cdi_singh')
col_order = sorted(alphasdf['study'].unique())
#col_order = [i if i != 'cdi_singh' else 'edd_singh' for i in col_order]

# ob_zhu is duplicate of nash_zhu, remove it.
col_order.remove('nash_zhu')
# get rid of the annoying OB-NASH label
alphasdf = alphasdf.replace('OB-NASH', 'NASH')

# re-order samples so that healthy controls will be plotted first
controls, cases = get_classes(alphasdf)
cases += ['nonCDI', 'postFMT_CDI']
alphasdf = pd.concat(
    [alphasdf.query('DiseaseState == @controls'),
     alphasdf.query('DiseaseState == @cases')])

labeldict = fmt.get_labeldict(col_order)
labeldict['ob_zhu'] = 'Zhu 2013, OB/NASH'

for g, subdf in alphasdf.groupby('alpha_metric'):
    fig = plot_alphadf(subdf, col_order, labeldict, g)
    plt.savefig(args.out + '.' + g + '.pdf')
