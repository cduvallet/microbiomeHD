#!/usr/bin/env python
"""
This script plots the AUCs from single-dataset classifiers vs.
those from healthy vs. disease classifiers.
"""

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Add util to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
import Formatting as fmt

def plot_aucs(aucs, x_col, y_col, groupby_col, colors):
    """
    Scatter plot aucs[x_col] vs aucs[y_col], colored by colors[groupby_col]

    Parameters
    ----------
    aucs : pandas DataFrame
        has x_col, y_col, and groupby_col
    x_col, y_col, groupby_col : str
    colors : dict
        values in groupby_col: color to plot
    """
    sns.set_style('white')
    fig, ax = plt.subplots(figsize=(4,3))
    ax.plot([0, 1], [0, 1], '--', c='0.95')
    ax.plot([0.5, 0.5], [0, 1], '--', c='0.95')
    ax.plot([0, 1], [0.5, 0.5], '--', c='0.95')
    for g, subdf in aucs.groupby(groupby_col):
        ax.scatter(subdf[x_col], subdf[y_col], c=colors[g], label=g.upper())
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    fig.tight_layout()
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    return fig, ax, lgd

p = argparse.ArgumentParser()
p.add_argument('single_dataset_rf', help='path to rf_results from '
    + ' single-dataset classifiers.')
p.add_argument('h_v_dis_rf', help='path to rf_results from general healthy '
    + 'vs. disease classifiers.')
p.add_argument('dataset_out_fig', help='file name to save leave-one-dataset-'
    + 'out figure to.')
p.add_argument('disease_out_fig', help='file name to save leave-one-disease-'
    + 'out figure to.')
args = p.parse_args()

# Get disease colors for plotting
disease_colors = fmt.get_disease_colors()

## Single-dataset results
fn_rf_results = args.single_dataset_rf
rf_results = pd.read_csv(fn_rf_results, sep='\t')
rf_results = rf_results.replace('edd_singh', 'cdi_singh')
rf_results['disease'] = rf_results['dataset'].map(lambda x: x.split('_')[0])

## Leave one dataset out
dataset_out = pd.read_csv(args.h_v_dis_rf, sep='\t')\
    .query('classifier == "dataset_out"')
aucs = pd.merge(rf_results[['dataset', 'roc_auc', 'disease']].drop_duplicates(),
                dataset_out[['dataset', 'auc', 'disease']].drop_duplicates())\
         .dropna()

fig, ax, lgd = plot_aucs(aucs, 'roc_auc', 'auc', 'disease', disease_colors)
ax.set_xlabel('AUC, single-dataset classifier')
ax.set_ylabel('AUC, healthy vs. disease classifier')
fig.savefig(args.dataset_out_fig, bbox_extra_artists=(lgd,),
    bbox_inches='tight')

## Leave one disease out
disease_out = pd.read_csv(args.h_v_dis_rf, sep='\t')\
    .query('classifier == "disease_out"')

# Average the original dataset AUCs by disease
mean_auc = pd.DataFrame(rf_results.groupby('disease').mean()['roc_auc'])
mean_auc = mean_auc.rename(columns={'roc_auc': 'mean_dataset_auc'})

# Concatenate averaged results with leave-disease-out results
aucs = disease_out[['disease', 'auc']].drop_duplicates()
aucs.index = aucs['disease']
aucs = aucs.rename(columns={'auc': 'auc_diseaseout'})
aucs = pd.concat((aucs, mean_auc), axis=1).dropna()
# Plot
fig, ax, lgd = plot_aucs(aucs, 'mean_dataset_auc', 'auc_diseaseout',
    'disease', disease_colors)
ax.set_xlabel('mean AUC, single-dataset classifier')
ax.set_ylabel('')
fig.savefig(args.disease_out_fig, bbox_extra_artists=(lgd,),
    bbox_inches='tight')
