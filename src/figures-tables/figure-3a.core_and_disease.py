#!/usr/bin/env python
"""
This script plots three panels: one with the disease-specific consistent
bugs, one with the core bugs, and one with the phylogeny.
"""
import argparse
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap # for phylo colors



def get_phylo_colors(keep_rows):
    """
    Return the df with phylum and orders. Make color_dict to plot bar phylogeny.

    Parameters
    ----------
    keep_rows         list of OTUs to keep. Should be full name (starting with k__;...)
                      and in the order that you'll want to plot them.

    Return
    ------
    phylodf           dataframe with 'phylum', 'class', 'order', and other columns.
                      values are text ('o__Clostridiales')
    palette           a ListedColormap of all the colors to map to phylum/orders
    color_dict        dictionary with {tax_level: index_in_palette}, where tax_level
                      is a string like 'o__Clostridiales'.
    """
    phylodf = pd.DataFrame(columns=['phylum', 'class', 'order', 'full'])
    phylodf['full'] = keep_rows
    phylodf['phylum'] = phylodf['full'].apply(lambda x: x.split(';')[1])
    phylodf['class'] = phylodf['full'].apply(lambda x: x.split(';')[2])
    phylodf['order'] = phylodf['full'].apply(lambda x: x.split(';')[3])
    phylodf['family'] = phylodf['full'].apply(lambda x: x.split(';')[4])
    phylodf['genus'] = phylodf['full'].apply(lambda x: x.split(';')[5])

    ## Set1 color palette
    colors = sns.color_palette('Set1', 9)
    # red (proteo), blue (bacteroides), green (actino), purple (firmicutes),
    # orange (fuso), yellow, brown (eury), pink (verruco), gray (teneri, cyano, lenti, synerg)
    # this color dict is for 'Set1' color palette
    color_dict = {'p__Actinobacteria': 2,
                  'p__Bacteroidetes': 1,
                  'p__Cyanobacteria/Chloroplast': 7,
                  'p__Candidatus_Saccharibacteria': 7,
                  'p__Euryarchaeota': 6,
                  'p__Firmicutes': 3,
                  'p__Fusobacteria': 4,
                  'p__Lentisphaerae': 7,
                  'p__Proteobacteria': 0,
                  'p__Synergistetes': 7,
                  'p__Verrucomicrobia': 8,
                  'p__Tenericutes': 7}

    ## Check that all phyla in phylodf have a color
    missingphyla = [i for i in phylodf['phylum'].unique() if i not in color_dict]
    if len(missingphyla) > 0:
        print('You need to give the following phyla colors in get_phylo_colors():')
        print('\n'.join(missingphyla))
    ## I think I also need to drop any phyla from the dict, otherwise my
    ## ListedColormap returns random-ish mappings? (It's not random, but I don't
    ## know how they map ...)
    drop_phyla = [i for i in color_dict if i not in phylodf['phylum'].unique()]
    for i in drop_phyla:
        color_dict.pop(i)
    # Seed a light palette from the top-level phylum color.
    # The first element in light_palette is white and the last element is the original color
    # So make size of palette = number of orders + 2, and pick the 2nd element to the 2nd to last
    order_color_dict = {}
    for p, subdf in phylodf.groupby('phylum'):
        # Get unique orders in that phylum
        orders = subdf['order'].unique()
        # Seed a light palette from the top-level phylum color.
        # The first element in light_palette is white and the last element is the original color
        # So make size of palette = number of orders + 2, and pick the 2nd element to the 2nd to last
        order_color_dict[p] = sns.light_palette(colors[color_dict[p]], len(orders) + 2)[1:-1]
        # Assign corresponding index in `colors` to each order, and add to the master color_dict
        color_dict.update({orders[i]: len(colors) + i for i in range(len(orders))})
        # Append the corresponding colors to the master list of colors
        colors += [tuple(order_color_dict[p][i]) for i in range(len(orders))]

    palette = ListedColormap(colors)

    return phylodf, palette, color_dict

p = argparse.ArgumentParser()
p.add_argument('disease', help='path to file with disease-wise meta-analysis '
               + 'results. Rows and columns should be ordered as you want '
               + 'them plotted.')
p.add_argument('overall', help='path to file with core/overall meta-analysis '
               + 'results. Rows and columns should be ordered as you want '
               + 'them plotted.')
p.add_argument('out', help='out file to save figure to')
p.add_argument('--labels', help='flag to add genus labels to plot. [default: '
               + '%(default)s]', action='store_true')
args = p.parse_args()

# Read in files
disease_meta = pd.read_csv(args.disease, sep='\t', index_col=0)
overall_meta = pd.read_csv(args.overall, sep='\t', index_col=0)

# Check that the rows are in the same order for both panels,
# otherwise raise error
if sum(overall_meta.index == disease_meta.index) != overall_meta.shape[0]:
    raise ValueError("Your overall and disease-wise meta-analysis "
                     + " dataframes don't contain the same rows in "
                     + "the same order!!")

# Prepare for plotting
phylodf, palette, color_dict = get_phylo_colors(disease_meta.index)
phylo_toplot = phylodf['order'].apply(lambda x: color_dict[x])

## Set up plot
sns.set_style('white')
# Set up 2 grid specs, one with the phylogeny and one with the two heatmaps
fig = plt.figure(figsize=(5, 16))

# Left most bar has order level heatmap
gsL = gridspec.GridSpec(nrows=1, ncols=1, left=0.45, right=0.51,
                        top=0.97, bottom=0.03)
# axL1 has the order-level phylogeny
axL1 = plt.subplot(gsL[0])

# Middle bar has two heatmaps (one with one column, one with four columns)
gsM = gridspec.GridSpec(nrows=1, ncols=5, left=0.56, wspace=0.6,
                        top=0.97, bottom=0.03)
# axM1 has the "overall" significant
axM1 = plt.subplot(gsM[0])
# axM2 has the disease-specific "overall" significant
axM2 = plt.subplot(gsM[1:])

## Phylogeny, left axis
# Manually get max and min values, otherwise colormap gets
# scaled to values within each column
vmax = phylo_toplot.max().max()
vmin = phylo_toplot.min().min()

axL1.imshow(pd.DataFrame(phylo_toplot), interpolation='nearest',
            cmap=palette, aspect='auto', vmax=vmax, vmin=vmin)
# Need to shift ylim a bit to not cutoff the edge cell
axL1.set_ylim(axL1.get_ylim()[0], axL1.get_ylim()[1] - 0.5)
axL1.set_xticklabels([])
[axL1.spines[i].set_visible(False) for i in
    ['right', 'top', 'left', 'bottom']]

## Meta-analysis, right axis
# values: 1 = disease = red, -1 = health = blue, 0 = mixed = black)
axM1.imshow(overall_meta.values, interpolation='nearest', aspect='auto',
            cmap=sns.diverging_palette(220,20,center='dark',as_cmap=True))
axM1.set_xticklabels([])
axM1.set_yticklabels([])
axM1.set_ylim(axM1.get_ylim()[0], axM1.get_ylim()[1] - 0.5)

axM2.imshow(disease_meta.values, interpolation='nearest', aspect='auto',
            cmap=sns.diverging_palette(220,20,center='dark',as_cmap=True))
axM2.set_yticklabels([])
axM2.set_xticklabels([])
axM2.set_ylim(axM2.get_ylim()[0], axM2.get_ylim()[1] - 0.5)

if args.labels:
    ## Label left-most axis with significant genera
    # Get series with True/False if any of the values are not NaN
    labeldf = overall_meta.join(disease_meta)
    labeldf = labeldf.applymap(np.isfinite).sum(axis=1).astype(bool)
    #labels = [i.split(';')[-1][3:] if labeldf.loc[i] else '' for i in labeldf.index]
    labels = [i.split(';')[-1][3:] for i in labeldf.index]

    #axL1.yaxis.tick_right()
    #axL1.yaxis.set_label_position("right")
    axL1.set_yticks(range(0, len(labeldf)))
    axL1.set_yticklabels(labels,fontsize='xx-small', rotation=180, va='bottom')

else:
    axL1.set_yticklabels([])

plt.savefig(args.out)
