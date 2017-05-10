#!/usr/bin/env python
"""
This script makes the supplementary figures showing results
from different random forest parameters.
"""
import argparse
import pandas as pd
import numpy as np

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
p.add_argument('rf_params', help='tidy file with RF parameter search results')
p.add_argument('criteria', help='criteria to plot', choices=['gini', 'entropy'])
p.add_argument('out', help='file to write figure to')
args = p.parse_args()

sns.set_style('white')

## Read in RF parameters
df = pd.read_csv(args.rf_params, sep='\t')
# cdi_youngster is nan, don't plot
df = df.dropna()

## Set up for plotting
df['est_leaf'] = df.apply(
    lambda row: str(row['n_estimators']) + '_' + str(row['min_samples_leaf']),
    axis=1)
hue_order = ['1000_1', '1000_2', '1000_3', '10000_1', '10000_2', '10000_3']
hue_kws = {'marker': ['v', 'v', 'v', '^', '^', '^']}
pal = {'1000_1': 'red', '10000_1': 'red',
       '1000_2': 'blue', '10000_2': 'blue',
       '1000_3': 'green', '10000_3': 'green'}

## Get order
dataset_order = df['dataset'].unique()
dataset_order = [i if i != 'edd_singh' else 'cdi_singh' for i in dataset_order]
dataset_order = np.sort(dataset_order)
dataset_order = [i if i != 'cdi_singh' else 'edd_singh' for i in dataset_order]

crit = args.criteria
g = sns.FacetGrid(data=df.query('criterion == @crit'), col='dataset',
                  col_wrap=6, hue='est_leaf', sharey=False, #ylim=(0,1),
                  palette=pal, hue_order=hue_order, hue_kws=hue_kws,
                  col_order=dataset_order, aspect=1.5, size=2)

g.map(sns.stripplot, 'min_samples_split', 'roc_auc', marker='D',
      size=10, jitter=True, alpha=0.75, edgecolor='gray')

labels = get_labeldict(dataset_order)
ax_counter = 0
for ax in g.axes:
    oldtitle = ax.get_title()
    newtitle = oldtitle.split('=')[1].strip()
    ax.set_title(labels[newtitle], fontsize='large')
    if ax_counter in [0, 6, 12, 18, 24]:
        ax.set_ylabel('AUC', fontsize='large')
    if ax_counter in range(23, 29):
        ax.set_xlabel('min_samples_split', fontsize='large')
    ax_counter += 1

fig = plt.gcf()
fig.set_tight_layout(True)
fig.savefig(args.out)
