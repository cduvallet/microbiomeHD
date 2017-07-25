#!/usr/bin/env python
"""
This script makes the disease-specific heatmaps of q-values, where each
genus (row) is sorted by overall disease-ness and health-ness.
"""
import argparse
import pandas as pd
import numpy as np
from string import upper

import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import seaborn as sns

def plot_disease_heatmap(disdf, samplesizes, vmax=None, with_labels=False):
    """
    Plot heatmap for one disease.

    Parameters
    ----------
    disdf : pandas dataframe
        contains signed q-values with ordered genera in rows and datasets
        in columns.
    samplesizes : pandas dataframe
        datasets in rows or 'dataset' column,
        and 'N_ctrl', 'N_dis', and 'total' columns
    vmax : float or int
        colorbar cuts off at +/- vmax (vmax should be positive)
    with_labels : bool (default: False)
        whether to plot labels or just return heatmap alone

    Returns
    -------
    fig : matplotlib Figure instance
    """

    if 'dataset' in samplesizes.columns:
        samplesizes.index = samplesizes['dataset']

    fig = plt.figure(figsize=(10,0.15*disdf.shape[0]))
    sns.set_style('white')
    ax_heatmap = plt.gca()

    ## Plot heatmap
    # Set vmin == -vmax
    if vmax is None:
        vmax = max(disdf.max().max(), -disdf.min().min())

    cax = ax_heatmap.imshow(disdf.values, interpolation='nearest',
                            aspect='auto',
                            cmap=sns.diverging_palette(220,20,as_cmap=True),
                            vmax=vmax, vmin=-vmax)
    ax_heatmap.set_aspect(0.1)

    if with_labels:
        ## Add dataset labels
        labels = []
        for dataset in disdf:
            label = dataset.split('_')[1]
            label = upper(label[0]) + label[1:]
            # Definitely a better way to do this but oh well...
            if dataset == 'cdi_schubert2':
                label = 'Schubert\n(nonCDI)'
            elif dataset == 'cdi_schubert':
                label = 'Schubert\n(CDI)'
            elif dataset in ['edd_singh', 'cdi_singh']:
                label = 'Singh\n(EDD)'
            elif dataset in ['cdi_vincent', 'cdi_youngster']:
                label += '\n(CDI)'
            label += '\n' + str(samplesizes.loc[dataset, 'N_ctrl']) \
                  + '\n' + str(samplesizes.loc[dataset, 'N_dis'])
            labels.append(label)
        ax_heatmap.set_xticks(range(disdf.shape[1]))
        ax_heatmap.set_xticklabels(labels, fontsize='large')

        # Make genus labels
        ax_heatmap.set_yticks(range(disdf.shape[0]))
        ax_heatmap.set_yticklabels([';'.join(i.split(';')[-2:]) for i in disdf.index], fontsize='x-small')
    else:
        ax_heatmap.set_yticks([])
        ax_heatmap.set_yticklabels([])
        ax_heatmap.set_xticks(range(0, disdf.shape[1]))
        ax_heatmap.set_xticklabels('')

    return fig

p = argparse.ArgumentParser()
p.add_argument('disease', help='disease to plot. Given disease should match '
               + 'dataset prefixes (e.g. "cdi", "hiv", etc)')
p.add_argument('qvalues', help='path to file with qvalues. Genera are in '
               + 'rows, and should be ordered as you want them plotted. '
               + 'Datasets are in columns, and will be ordered by sample '
               + 'size.')
p.add_argument('dataset_info', help='path to file with info about datasets. '
               + 'Should have columns "dataset", "H", "case", and "total".')
p.add_argument('out_file', help='path to write figure file to')
p.add_argument('--labels', help='flag to include genus and dataset labels '
               + 'on plot. [default: %(default)s]', action='store_true')
p.add_argument('--qthresh', help='q-value that should be the most opaque '
               + 'in the heatmap. Also used as the significance threshold '
               + 'when keeping only rows with at least one significant genus. ' + '[default: %(default)s]', default=0.05)
args = p.parse_args()

# Read in pvalues. df has datasets in columns, genera in rows
df = pd.read_csv(args.qvalues, sep='\t', index_col=0)
# Replace edd_singh with cdi_singh (for text parsing reasons)
df = df.rename(columns={'edd_singh': 'cdi_singh',
                        'noncdi_schubert': 'cdi_schubert2'})
# Name the index 'otu' for future melting etc
df.index.name = 'otu'

samplesizes = pd.read_csv(args.dataset_info, sep='\t')
samplesizes = samplesizes\
    .replace('edd_singh', 'cdi_singh')\
    .replace('noncdi_schubert', 'cdi_schubert2')
samplesizes = samplesizes.sort_values(by='total', ascending=False)

## Replace pvalues with +/- 1 if significant or not
sigmap = lambda x: np.sign(x) if abs(x) < args.qthresh else 0
dfsig = df.applymap(sigmap)
datasets = df.columns

disease = args.disease + '_'
keep_datasets = [i for i in datasets if i.startswith(disease)]

disdf = dfsig[keep_datasets]
# Keep only OTUs which were signficant in at least one study
disdf = df.loc[disdf.apply(abs).sum(axis=1) != 0, keep_datasets]

# Convert to signed log10(pval)
disdf = disdf.applymap(lambda x: np.sign(x)*abs(np.log10(abs(x))))

# Re-order rows
disdf = disdf.loc[disdf.sum(axis=1).sort_values(ascending=False).index]
# Re-order columns by sample size
disdf = disdf[[i for i in samplesizes['dataset'] if i in keep_datasets]]
fig = plot_disease_heatmap(disdf, samplesizes, vmax=abs(np.log10(args.qthresh)),
                           with_labels=args.labels)
fig.tight_layout()
fig.savefig(args.out_file)
