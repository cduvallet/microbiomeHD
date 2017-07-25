#!/usr/bin/env python
"""
This script makes Figure 1, which has the sample sizes, AUC, number of
significant genera, and direction of shift.
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
import Formatting as fmt


def stemplot(x, y, data, order, ax, palette, marker='o', size=7):
    """
    Wrapper to make one stemplot with colored dashed lines leading
    to colored marker.

    Parameters
    ----------
    x, y : str
        used in call to sns.stripplot() with data
    data : pandas dataframe
        Should have the values that are given for 'order' in the index,
        or a column called 'label' with those values.
    order : list
        order of x values
    ax : Axis object
        axis handle to plot values on
    palette : dict
        {values in x-axis : color mapping value}
    marker : str
        marker value to pass to stripplot
    size : int
        size of marker

    Returns
    -------
    ax
    """

    if 'label' in data:
        data.index = data['label']

    sns.stripplot(x=x, y=y, data=data, order=order, ax=ax, palette=palette,
                  size=size, marker=marker)
    _, stemlines, baseline = ax.stem(data.loc[order, y],
                                     markerfmt=" ", linefmt=":")
    # Remove stemplot baseline
    plt.setp(baseline, visible=False)

    # Change stem colors
    colorslist = [palette[i] for i in order]
    _ = [plt.setp(stemlines[i], 'color', colorslist[i])
            for i in range(len(colorslist))]
    _ = [i.set_alpha(0.75) for i in stemlines]

    return ax

def plot_fig1(dysbiosis, dataset_order, samplesizes, edd_color=False):
    """
    Make figure 1 with sample sizes, AUCs, number of significant genera,
    and direction of shift.

    Note that 'edd_singh' should be converted to 'cdi_singh' in all of
    these inputs.

    Parameters
    ----------
    dysbiosis : pandas DataFrame
        should have 'metric', 'label', 'value' columns, where
        'metric' includes at least ['auc', 'n_sig', 'balance'] and 'label'
        column contains the dataset labels.
    dataset_order : list
        order of datasets to plot on the x-axis.
    samplesizes : pandas dataframe
        df with at least 'dataset' and 'total' columns
    edd_color : bool
        whether to separate the edd study with a manually defined color
    """
    disease_colors = fmt.get_disease_colors()
    # Make color palette dictionary that has all of the datasets
    # Note: need trailing underscore so that cd studies are not considered cdi
    diseases = set([i.split('_')[0] + '_' for i in dataset_order])
    colors = {}
    for d in diseases:
        dis_datasets = [i for i in dataset_order if i.startswith(d)]
        colors.update({i: j for i, j in zip(dis_datasets,
            len(dis_datasets)*[sns.light_palette(disease_colors[d[:-1]])[-1]])})
        if edd_color:
            colors['cdi_singh'] = disease_colors['edd']

    ### Plot
    sns.set_style('white', {'ytick.direction': 'in', 'ytick.major.size': 2.0,
                           'font.size': 'xx-small'})

    ## Set up the plot
    # For spacing between plots: hspace in the following calls to gridspec
    # set the space between the subplots.
    # To set the space between the two gridspecs (top and bottom), need to
    # change h_pad in call to tight_layout() at the very bottom!
    fig = plt.figure(figsize=(0.27*len(dataset_order),8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    # Bottom gridspec has sample size and AUCs
    gs_bottom = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1],
                                                 hspace=0.25)
    ax_samplesize = plt.subplot(gs_bottom[1])
    ax_auc = plt.subplot(gs_bottom[0])

    # Top gridspec has n_sig and balance
    gs_top = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0],
                                              hspace=0.25)
    ax_metrics = [plt.subplot(i) for i in gs_top]

    ### Plot stuff
    ## Sample size
    # Sample size is the one with the labels.
    # dataset_order = [i if i != 'cdi_singh' else 'edd_singh' for i in dataset_order]
    # samplesizes['dataset'] = [i if i != 'cdi_singh' else 'edd_singh' for i in samplesizes['dataset']]
    # colors = {(i if i != 'cdi_singh' else 'edd_singh'): colors[i] for i in colors}
    labeldict = fmt.get_labeldict(dataset_order)

    # Plot barplot.
    sns.barplot(x='dataset', y='total', data=samplesizes,
                order=dataset_order, ax=ax_samplesize, palette=colors,
                edgecolor='')

    [ax_samplesize.spines[i].set_visible(False) for i in ['left', 'top']]

    ax_samplesize.yaxis.tick_right()
    ax_samplesize.yaxis.set_label_position('right')
    ax_samplesize.set_yticks([100,300,500])
    ax_samplesize.set_yticklabels(['100', '300', '500'], rotation=90)

    ax_samplesize.set_xticklabels([labeldict[i] for i in dataset_order],
                                  rotation=90)

    ax_samplesize.set_xlabel('')
    ax_samplesize.set_ylabel('Total samples')

    # # Change all labels back to the original cdi_singh
    # dataset_order = [i if i != 'edd_singh' else 'cdi_singh' for i in dataset_order]
    # samplesizes['dataset'] = [i if i != 'edd_singh' else 'cdi_singh' for i in samplesizes['dataset']]
    # colors = {(i if i != 'edd_singh' else 'cdi_singh'): colors[i] for i in colors}

    ## AUCs
    ## Plot AUCs as barplots
    # Plot dashed line
    ax_auc.plot([-0.5, 29.5], [1, 1], color='0.75', linestyle='--')

    # No dashed lines
    # sns.stripplot(x='label', y='value',
    #               data=dysbiosis.query('metric == "auc"'),
    #               order=dataset_order, ax=ax_auc, palette=colors,
    #               size=7, marker='D')
    # I liked it with the dashed lines better
    ax_auc = stemplot(x='label', y='value',
                      data=dysbiosis.query('metric == "auc"'),
                      order=dataset_order, ax=ax_auc,
                      palette=colors, marker='D')

    ax_auc.set_xticklabels([])
    ax_auc.set_ylim([0.5, 1])
    ax_auc.set_yticks([0.5, 0.75, 1.0])
    ax_auc.set_yticklabels([0.5, 0.75, 1], rotation=90)
    ax_auc.yaxis.tick_right()
    ax_auc.yaxis.set_label_position('right')

    [ax_auc.spines[i].set_visible(False) for i in ['left', 'top']]

    ax_auc.set_ylabel('AUC')
    ax_auc.set_xlabel('')

    ## n_sig
    ax_metrics[1] = stemplot(x='label', y='value',
                             data=dysbiosis.query('metric == "n_sig"').replace(0, np.nan),
                             order=dataset_order, ax=ax_metrics[1],
                             palette=colors)#, marker='o')

    ax_metrics[1].set_xticklabels([])
    ax_metrics[1].set_ylim([0, None])
    ax_metrics[1].set_yticks([0, 20, 40, 60])
    ax_metrics[1].set_yticklabels([0, 20, 40, 60], rotation=90)
    ax_metrics[1].yaxis.tick_right()
    ax_metrics[1].yaxis.set_label_position('right')

    [ax_metrics[1].spines[i].set_visible(False) for i in ['left', 'top']]

    ax_metrics[1].set_ylabel('Genera, q < 0.05')
    ax_metrics[1].set_xlabel('')

    ## balance
    # center the balance values around zero
    dysbiosis['centered_value'] = \
        dysbiosis.query('metric == "balance"')['value'] - 0.5

    ax_metrics[0] = stemplot(x='label', y='centered_value',
                             data=dysbiosis.query('metric == "balance"'),
                             order=dataset_order, ax=ax_metrics[0],
                             palette=colors)

    # Plot the gray, blue and red lines for 0% and 100% depleted/enriched
    ax_metrics[0].plot([-0.5, 29.5], [0, 0],
                       color='0.75', linestyle='-')
    ax_metrics[0].plot([-0.5, 29.5], [-0.5, -0.5],
                       color='b', linestyle='-', alpha=0.2)
    ax_metrics[0].plot([-0.5, 29.5], [0.5, 0.5],
                       color='r', linestyle='-', alpha=0.2)

    ax_metrics[0].set_ylim([-0.6, 0.6])
    ax_metrics[0].yaxis.tick_right()
    ax_metrics[0].yaxis.set_label_position('right')
    # The -0.5 tick needs to be moved a bit to be less scrunched
    ax_metrics[0].set_yticks([-0.4, 0.5])
    ax_metrics[0].tick_params(axis='y', which=u'both', length=0)
    ax_metrics[0].set_yticklabels(['Beneficial-\n depleted',
                                   'Pathogen-\n enriched'],
                                  rotation=90)
    ax_metrics[0].set_xticklabels([])
    ax_metrics[0].yaxis.set_label_position('right')

    [ax_metrics[0].spines[i].set_visible(False) for i in ['left', 'top', 'bottom']]

    ax_metrics[0].set_xlabel('')
    ax_metrics[0].set_ylabel('')

    # Tight layout and set spacing between A and B panels
    gs.tight_layout(fig, h_pad=3.5)

    return fig

p = argparse.ArgumentParser()
p.add_argument('dysbiosis', help='path to tidy table with dysbiosis metrics')
p.add_argument('dataset_info', help='path to table with samplesize info')
p.add_argument('out_file', help='file to save figure as')
p.add_argument('--edd', help='flag to color edd_singh differently than cdi',
    action='store_true')
args = p.parse_args()

dysbiosis = pd.read_csv(args.dysbiosis, sep='\t')
dataset_info = pd.read_csv(args.dataset_info, sep='\t')

# Note: need to replace edd_singh with cdi_singh in basically everything.
dysbiosis = dysbiosis\
    .replace('edd_singh', 'cdi_singh')\
    .replace('noncdi_schubert', 'cdi_schubert2')
dataset_info = dataset_info\
    .replace('edd_singh', 'cdi_singh')\
    .replace('noncdi_schubert', 'cdi_schubert2')

_, dataset_order = fmt.get_dataset_order(dataset_info)

fig = plot_fig1(dysbiosis, dataset_order, dataset_info, args.edd)
fig.savefig(args.out_file)
