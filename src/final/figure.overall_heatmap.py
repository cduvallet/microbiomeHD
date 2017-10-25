#!/usr/bin/env
"""
This script plots the values given for each genus in each dataset.
"""

import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('TKAgg')
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import argparse


# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
import Formatting as fmt


def plot_overall_heatmap_figure(mean_toplot, phylo_toplot, overall_df,
                                disease_df, dataset_order,
                                heatmap_vmax=None, figsize=(13,14)):
    """
    Heatmap with all datasets + overall results.
    Left panel: color bars corresponding to phylogeny
    Second panel: color bar with overall significant genera
    Third panel: heatmap with sig genera per disease
    Right panel: p-values (or effect) for each dataset

    No re-ordering of rows or columns is done.

    Parameters
    ----------
    mean_toplot : pandas DataFrame
        df with values to plot in right panel. Genera (rows) should already be
        ordered as they will be plotted. Columns are re-ordered by the order in
        dataset_order, and grouped by disease (i.e. prefix in each dataset in
        dataset_order). edd_singh should be cdi_singh in this dataframe.
    phylo_toplot : numpy array
        Values to plot for phylogeny heatmap. Values should be tuples,
        corresponding RGBA values.
    overall_df : pandas DataFrame or Series
        df with one column, with 1 (disease), -1 (health), and 0 (both)
        values. Values will be plotted red, blue, and black (respectively).
    disease_df : pandas DataFrame
        dataframe with one column per disease type, and 1/-1/0 values,
        as in overall_df
    dataset_order : list
        columns in disease_df in the order that they should be plotted.
        Note that prefixes should be grouped by disease (the datasets
        are split on '_' delimiter and consecutive datasets with the same
        prefix are plotted in the same sub-heatmap)
    heatmap_vmax : float
        value to use as vmax for the dataset heatmap
    """
    ## Default gridspec axis locations are here: http://matplotlib.org/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels
    #left = 0.125    #the left side of the subplots of the figure
    #right = 0.9    #the right side of the subplots of the figure
    #bottom = 0.1    #the bottom of the subplots of the figure
    #top = 0.9    #the top of the subplots of the figure
    #wspace = 0.2    #the amount of width reserved for blank space between subplots
    #hspace = 0.2    #the amount of height reserved for white space between subplots

    ## Set up 3 grid specs on my plot
    fig = plt.figure(figsize=figsize)
    # Left bar has one heatmap (order level)
    gsL = gridspec.GridSpec(1, 1)
    gsL.update(left=0.12, right=0.14, bottom=0.15, top=0.98)
    # Middle bar has two heatmaps (one with one column, one with four columns)
    gsM = gridspec.GridSpec(1,5)
    gsM.update(left=0.16, right=0.27, wspace=0.6, bottom=0.15, top=0.98)
    # Right bar has a bunch of heatmaps (one for each disease)
    gsR = gridspec.GridSpec(1, mean_toplot.shape[1])
    gsR.update(left=0.28, right=0.86, wspace=0.3, bottom=0.15, top=0.98)

    ## Phylogeny, left axis
    # axL1 has order-level colors (only column in phylodf)
    axL1 = plt.subplot(gsL[0])
    axL1.imshow(phylo_toplot, interpolation='nearest', aspect='auto')
    axL1.set_axis_off()
    axL1.set_ylim(axL1.get_ylim()[0], axL1.get_ylim()[1] - 0.5)

    ## Meta-analysis, middle axis
    # axM1 has the "overall" significant
    # values: 1 = disease = red, -1 = health = blue, 0 = mixed = purple)
    axM1 = plt.subplot(gsM[0])
    axM1.imshow(overall_df.values, interpolation='nearest', aspect='auto',
        cmap=sns.diverging_palette(220,20,center='dark',as_cmap=True))
    #axM1.set_xticklabels([])
    axM1.set_xticks([0.4])
    axM1.set_xticklabels(['Non-specific'], rotation=45,
        fontsize='small', ha='right')
    axM1.set_yticklabels([])
    axM1.set_ylim(axM1.get_ylim()[0], axM1.get_ylim()[1] - 0.5)
    # axM2 has the disease-specific "overall" significant
    axM2 = plt.subplot(gsM[1:])
    axM2.imshow(disease_df.values, interpolation='nearest', aspect='auto',
        cmap=sns.diverging_palette(220,20,center='dark',as_cmap=True))
    axM2.set_yticklabels([])
    axM2.set_xticks(np.arange(0.4, disease_df.shape[1] + 0.4))
    labels = disease_df.rename(columns={'CDI': 'Diarrhea'})
    axM2.set_xticklabels(labels, rotation=45,
        fontsize='small', ha='right')
    axM2.set_ylim(axM2.get_ylim()[0], axM2.get_ylim()[1] - 0.5)

    ## Dataset-specific values, right axis
    # Get the total number of diseases, and the indices of each
    alldis = []
    for d in dataset_order:
        d = d.split('_')[0]
        if d not in alldis:
            alldis.append(d)

    if heatmap_vmax is None:
        heatmap_vmax = max(mean_toplot.max().max(), -mean_toplot.min().min())

    for dis in alldis:
        indices = [i for i in range(len(dataset_order))
                   if dataset_order[i].startswith(dis)]
        # Size each axis (i.e. number of gridspec slots) according to the
        # number of studies in that disease
        ax = plt.subplot(gsR[:, indices[0]:indices[0]+len(indices)])
        subdf = mean_toplot[[i for i in mean_toplot.columns
                             if i.startswith(dis)]]
        _ = ax.imshow(subdf.values, interpolation='nearest', aspect='auto',
                      cmap=sns.diverging_palette(220,20,as_cmap=True),
                      vmax=heatmap_vmax, vmin=-heatmap_vmax)
        ax.set_ylim([ax.get_ylim()[0], ax.get_ylim()[1] - 0.5])
        ax.set_yticklabels([])
        ax.set_xticks(range(0, subdf.shape[1]))
        # Align rotated labels with actual bar
        ax.set_xticks(ax.get_xticks() + 0.4)
        labels = fmt.get_labeldict(subdf.columns)
        labels = [labels[i] for i in subdf.columns]
        ax.set_xticklabels(labels, rotation=45, fontsize='small',
                           ha='right')
        ax.spines['bottom'].set_color('gray')
        ax.spines['left'].set_color('gray')
        ax.spines['right'].set_color('gray')
        ax.spines['top'].set_color('gray')

    # Add genus labels to right-most axis
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')
    labels = [i.split(';')[-1][3:] for i in subdf.index]
    ax.set_yticks(np.arange(0, len(labels)))
    ax.set_yticklabels(labels, fontsize=8, va='center')

    return fig

def prepare_heatmap_plot(main_values, disease_meta, overall_meta, dataset_info,
                         plot_log10=False, qthresh=0.05):
    """
    Prepare the dataframes and colors for plotting the large heatmap.

    Parameters
    ----------
    main_values : pandas DataFrame
        DataFrame with values to plot for each genus in each dataset.
        Should have datasets in columns in any order and genera in rows in the
        order to be plotted. Values should be signed according to effect
        direction.
    disease_meta : pandas DataFrame
        DataFrame with disease-specific meta analysis, exact values to plot
        (i.e. +/- 1, no transformation of these values is done before
        plotting). Rows should be in the same order as main_values.
    overall_meta : pandas DataFrame
        DataFrame with overall significant bugs, values to plot (i.e. 1, -1, 0)
        Rows should be in the same order as main_values.
    dataset_info : pandas DataFrame
        DataFrame with sample sizes for all datasets. Passed into function
        get_dataset_order().
    plot_log10 : True/False (default: False)
        Whether to transform the p-values in the datasetheatmap to
        log10(p-values) (matching the original p-value's sign)
    qthresh : float
        q-value threshold to use as the most opaque color. Only set if
        plot_log10 is True.

    Returns
    -------
    main_values, phylo_toplot, overall_meta, disease_meta,
    dataset_order, vmax
        Inputs to plot_overall_heatmap_figure()
    """

    ## Prepare dataset orders
    dataset_info = dataset_info\
        .replace('edd_singh', 'cdi_singh')\
        .replace('noncdi_schubert', 'cdi_schubert2')
    _, dataset_order = fmt.get_dataset_order(dataset_info)
    # Need to convert edd_singh to cdi_singh
    main_values = main_values.rename(columns={'edd_singh': 'cdi_singh',
        'noncdi_schubert': 'cdi_schubert2'})

    ## Prepare phylogeny colors
    phylodf, color_dict = fmt.get_phylo_colors(main_values.index)
    # Convert phylum and order columns to numbers corresponding to color
    # from the palette, which has RGBA values
    phylo_toplot = phylodf['order'].map(lambda x: color_dict[x])
    phylo_toplot = np.array([[i] for i in phylo_toplot.values])

    ## Plot
    # Re-order columns of main dataframe according to dataset_order
    if main_values.shape[1] != len(dataset_order):
        print("check hard-coding of disease_order in get_dataset_order() - " +
              "it looks like you're missing some datasets!")
    main_values = main_values[dataset_order]

    # If needed, transform and set vmax
    vmax = None
    if plot_log10:
        logmap = lambda x: np.sign(x)*abs(np.log10(abs(x)))
        main_values = main_values.applymap(logmap)
        vmax = abs(np.log10(qthresh))

    return (main_values, phylo_toplot, overall_meta, disease_meta,
           dataset_order, vmax)

p = argparse.ArgumentParser()
p.add_argument('main_values', help='file with values for each dataset to plot. '
    + 'Rows should be ordered as they will be plotted; columns can be in any '
    + 'order.')
p.add_argument('disease_meta', help='file with the disease-specific bacteria. '
    + 'Genera should be in rows, with diseases in the columns. Neither rows '
    + 'nor columns are re-ordered before plotting.')
p.add_argument('overall', help='file the the "core" genera. Genera are in '
    + 'rows, and file has one column labeled "overall". Rows are not '
    + 're-ordered before plotting.')
p.add_argument('dataset_info', help='file with samplesize info. Should have '
    + 'columns "total" and "dataset".')
p.add_argument('figout', help='file name to save figure to')
p.add_argument('--plot-log10', help='Flag to plot log10(values). [default: '
    + '%(default)s]', action='store_true')
p.add_argument('--qthresh', help='q-value threshold to use as most opaque '
    + 'in big heatmap. [default: %(default)s]', default=0.05, type=float)
p.add_argument('--figsize', help='comma-separated values for figure size. '
    + ' e.g. should be entered as "--figsize 13,14" for figure with width 13 '
    + ' height 14', default='13,14')
args = p.parse_args()


# Read in files
main_values = pd.read_csv(args.main_values, sep='\t', index_col=0)
disease_meta = pd.read_csv(args.disease_meta, sep='\t', index_col=0)
overall_meta = pd.read_csv(args.overall, sep='\t', index_col=0)
dataset_info = pd.read_csv(args.dataset_info, sep='\t')
figsize = tuple(int(i) for i in args.figsize.split(','))

if not (main_values.index == disease_meta.index).all() and \
       (main_values.index == overall_meta.index).all():
    raise ValueError("Rows in the three heatmaps are not the same!")

# capitalize disease abbreviations
disease_meta.columns = [i.upper() for i in disease_meta.columns]

for_plotting = prepare_heatmap_plot(main_values, disease_meta,
                                    overall_meta, dataset_info,
                                    plot_log10=args.plot_log10,
                                    qthresh=args.qthresh)

# Plot!
sns.set_style('white')
fig = plot_overall_heatmap_figure(*for_plotting, figsize=figsize)
fig.savefig(args.figout)
