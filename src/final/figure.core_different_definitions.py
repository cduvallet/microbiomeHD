#!/usr/bin/env python
"""
This script plots three types of core bugs: the normal heuristic (2 diseases),
the heuristic minus diarrhea datasets, and Stouffer's results.
"""
import argparse
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('TKAgg')
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
from Formatting import get_phylo_colors, reorder_index_from_tree

p = argparse.ArgumentParser()
p.add_argument('heuristic', help='path to file with core/overall meta-analysis '
    + 'results. Rows should be ordered as you want them plotted.')
p.add_argument('diarrhea', help='path to file with core bugs defined without '
    + 'diarrhea datasets.')
p.add_argument('stouffer', help='path to file with core bugs defined using '
    + 'Stouffers method.')
p.add_argument('tree_file', help='path to newick tree to order rows')
p.add_argument('out', help='out file to save figure to')
p.add_argument('--labels', help='flag to add genus labels to plot. [default: '
    + '%(default)s]', action='store_true')
args = p.parse_args()

# Read in files
cdi_core = pd.read_csv(args.heuristic, sep='\t', index_col=0)
nocdi_core = pd.read_csv(args.diarrhea, sep='\t', index_col=0)
stouffer_core = pd.read_csv(args.stouffer, sep='\t', index_col=0).dropna()
cdi_core.columns = ['2 diseases']
nocdi_core.columns = ['No diarrhea']
stouffer_core.columns = ['Stouffer']

# This re-orders rows alphabetically I think
core_all = pd.concat((stouffer_core, nocdi_core, cdi_core), axis=1)
# Reorder phylogenetically
ordered_rows = reorder_index_from_tree(args.tree_file, core_all.index)
core_all = core_all.loc[ordered_rows]

# Prepare for plotting
phylodf, color_dict = get_phylo_colors(core_all.index)
phylo_toplot = phylodf['order'].apply(lambda x: color_dict[x])
# Plot RGBA directly by making a n_genera X 1 X 4 array
phylo_toplot = np.array([[i] for i in phylo_toplot.values])

### Set up plot
sns.set_style('white')
fig = plt.figure(figsize=(2.5, 15))

# Left most bar has order level heatmap

# gsL = gridspec.GridSpec(1,1)
# gsL.update(left=0.12, right=0.35, bottom=0.15, top=0.95)
gsL = gridspec.GridSpec(nrows=1, ncols=1, left=0.45, right=0.51,
                        top=0.9, bottom=0.03)
axL1 = plt.subplot(gsL[0])
axL1.imshow(phylo_toplot, interpolation='nearest', aspect='auto')
axL1.set_ylim(axL1.get_ylim()[0], axL1.get_ylim()[1] - 0.5)
[axL1.spines[i].set_visible(False) for i in
    ['right', 'top', 'left', 'bottom']]
axL1.set_xticklabels([])

# gsM = gridspec.GridSpec(1,1)
# gsM.update(left=0.4, right=0.95, wspace=0.6, bottom=0.15, top=0.95)
gsM = gridspec.GridSpec(nrows=1, ncols=1, left=0.56, wspace=0.6,
                        top=0.9, bottom=0.03)

axM1 = plt.subplot(gsM[0])
axM1.imshow(core_all.values, interpolation='nearest', aspect='auto',
    cmap=sns.diverging_palette(220,20,center='dark',as_cmap=True))
axM1.set_yticklabels([])
axM1.set_ylim(axM1.get_ylim()[0], axM1.get_ylim()[1] - 0.5)

axM1.xaxis.tick_top()
axM1.set_xticks(range(0, core_all.shape[1]))
axM1.set_xticklabels(['Stouffer', 'No diarrhea', '2+ diseases'], rotation=270)

if args.labels:
    ## Label left-most axis with significant genera
    # Get series with True/False if any of the values are not NaN
    labels = [i.split(';')[-1][3:] for i in core_all.index]
    axL1.set_yticks(range(0, core_all.shape[0]))
    axL1.set_yticklabels(labels,fontsize='xx-small', rotation=180, va='center')

else:
    axL1.set_yticklabels([])

plt.savefig(args.out)
