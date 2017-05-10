#!/usr/bin/env python
"""
This script plots the percent overlap of each dataset's associations
with the "core" response.
"""
import argparse
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
import Formatting as fmt

p = argparse.ArgumentParser()
p.add_argument('dysbiosis', help='file with tidy dysbiosis metrics, including '
               + 'at least "perc_overlap".')
p.add_argument('dataset_info', help='file with info about datasets. Should '
               + 'columns "total" and "dataset".')
p.add_argument('out', help='file to write percent overlap figure to.')
args = p.parse_args()

# Prepare data for plotting
dysbiosis = pd.read_csv(args.dysbiosis, sep='\t')
samplesizes = pd.read_csv(args.dataset_info, sep='\t')
_, dataset_order = fmt.get_dataset_order(samplesizes)

toplot = dysbiosis.query('metric == "perc_overlap"')
toplot['value'] = toplot['value'].astype(float)
toplot = toplot.dropna()
neworder = [i for i in dataset_order if i in toplot['label'].tolist()]

labeldict = fmt.get_labeldict_for_overlap(neworder)

# Plot
sns.set_style('white', {'ytick.direction': 'in', 'ytick.major.size': 0.0})
fig, ax = plt.subplots(figsize=(6,3))

sns.barplot(x='label', y='value', data=toplot, order=neworder, ax=ax,
            color="0.65")

# Align rotated labels with actual bar
ax.set_xticks(ax.get_xticks() + 0.4)
ax.set_xticklabels([labeldict[i] for i in neworder], rotation=45,
                   ha='right')

ax.set_ylabel("Fraction overlap with core")
ax.set_yticks([0,0.5,1])
ax.set_xlabel('')
[ax.spines[i].set_visible(False) for i in ['right', 'top']]
fig.tight_layout()

fig.savefig(args.out)
