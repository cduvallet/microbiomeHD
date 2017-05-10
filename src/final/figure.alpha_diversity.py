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

p = argparse.ArgumentParser()
p.add_argument('alpha', help='file with tidy alpha diversity data')
p.add_argument('out', help='path to write alpha diversity figure to')
args = p.parse_args()

## Prepare alpha diversity dataframe
alphasdf = pd.read_csv(args.alpha, sep='\t')

# Remove un-annotated samples
alphasdf = alphasdf.replace(' ', np.nan)
alphasdf = alphasdf.replace('IBDundef', np.nan)
alphasdf = alphasdf.dropna()

col_order = sorted([i if i != 'edd_singh' else 'cdi_singh'
                    for i in alphasdf['study'].unique()])
col_order = [i if i != 'cdi_singh' else 'edd_singh' for i in col_order]
# ob_zhu is duplicate of nash_zhu
col_order.remove('nash_zhu')
# get rid of the annoying OB-NASH label
alphasdf = alphasdf.replace('OB-NASH', 'NASH')

labeldict = fmt.get_labeldict(col_order)
labeldict['ob_zhu'] = 'Zhu 2013, OB/NASH'

sns.set_style('white')
g = sns.FacetGrid(alphasdf, col='study', col_wrap=6,
                  col_order=col_order, sharex=False, sharey=False)
g = g.map(sns.boxplot, "DiseaseState", "SDI")
g = g.map(sns.stripplot, "DiseaseState", "SDI", split=True, jitter=True,
          size=5, linewidth=0.6)

fig = plt.gcf()
fig.set_size_inches(15.5, 8.75)

# Fix y-axis gridlines
for ax in g.axes:
    yticks = ax.get_yticks()
    # Get rid of 0.5 values...
    ax.set_ylim(floor(yticks[0]), ceil(yticks[-1]))
    if yticks[0] < 0:
        ax.set_ylim(0, None)
    yticks = ax.get_yticks()
    ax.set_yticks(yticks[::2])
    # Need some space on the y-axis for p-values
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]+1)
    # Update title
    oldtitle = ax.get_title()
    newtitle = labeldict[oldtitle.split('=')[1].strip()]
    ax.set_title(newtitle)

plt.tight_layout()
plt.savefig(args.out)
