#!/usr/bin/env python
"""
This script plots the ROC curves for all the datasets.
"""
import argparse

import pandas as pd

import matplotlib
matplotlib.use('TKAgg')
import seaborn as sns
import matplotlib.pyplot as plt

# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
from Formatting import get_labeldict

p = argparse.ArgumentParser()
p.add_argument('rf_results', help='tidy RF results')
p.add_argument('out', help='file to write figure to')
args = p.parse_args()

## Plot Faceted ROC curves
rfresults = pd.read_csv(args.rf_results, sep='\t')
rfresults = rfresults.query('dataset != "cdi_youngster"')

# Sort datasets alphabetically but group edd_singh and noncdi_schubert
# with the CDIs
col_order = sorted(
    rfresults['dataset']\
        .replace('edd_singh', 'cdi_singh')\
        .replace('noncdi_schubert', 'cdi_schubert2')\
    .unique())
col_order = [i if i != 'cdi_singh' else 'edd_singh' for i in col_order]
col_order = [i if i != 'cdi_schubert2' else 'noncdi_schubert'
             for i in col_order]

sns.set_style('white')
g = sns.FacetGrid(data=rfresults, col='dataset', col_wrap=6,
                  col_order=col_order, size=2, aspect=1)
g.map(plt.plot, 'mean_fpr', 'mean_tpr')

labels = get_labeldict(col_order)
ax_counter = 0
for ax in g.axes:
    ax.plot([0, 1], [0, 1], linestyle='--', color='0.75')
    oldtitle = ax.get_title()
    newtitle = oldtitle.split('=')[1].strip()
    ax.set_title(labels[newtitle])
    if ax_counter in [0, 6, 12, 18, 24]:
        ax.set_ylabel('TPR', fontsize='medium')
    if ax_counter in range(23, 29):
        ax.set_xlabel('FPR', fontsize='medium')
    ax_counter += 1

plt.savefig(args.out)
