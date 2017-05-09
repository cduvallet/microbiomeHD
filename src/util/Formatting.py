#!/usr/bin/env python
"""
Useful functions in formatting tables and figures.
"""
import numpy as np
import pandas as pd
import seaborn as sns

def get_dataset_order(df):
    """
    Given a list of diseases and a dataframe with 'total' sample size
    in a column, return the order of datasets.
    For each disease (in the order given), datasets are ordered from
    largest sample size to smallest.

    Parameters
    ----------
    df : pandas dataframe
        should have columns 'dataset' and 'total', where 'dataset' starts
        with the disease state ('edd' should be changed to 'cdi' for correct
        grouping)

    Returns
    -------
    disease_order : list
        hard-coded order of diseases
    dataset_order : list
        ordered values in 'dataset' column, ordered first by disease and
        then by descending 'total' sample size
    """
    df = df.sort_values(by='total', ascending=False)

    ## Set up dataset order
    disease_order = ['crc', 'ob', 'ibd', 'cdi', 'hiv', 'asd', 't1d', 'nash',
                     'mhe', 'ra', 'par']
    dataset_order = np.concatenate(
        [df.loc[
            df['dataset'].apply(lambda x:  x.startswith(d)), 'dataset'
               ].values
        for d in disease_order])

    return disease_order, dataset_order

def get_labeldict(dataset_order):
    d = {'asd_kang': 'Kang 2013, ASD',
         'asd_son': 'Son 2015, ASD',
         'cdi_schubert': 'Schubert 2014, CDI',
         'cdi_singh': 'Singh 2015, EDD',
         'cdi_vincent': 'Vincent 2013, CDI',
         'cdi_youngster': 'Youngster 2014, CDI',
         'crc_baxter': 'Baxter 2016, CRC',
         'crc_chen': 'Chen 2012, CRC',
         'crc_wang': 'Wang 2012, CRC',
         'crc_zackular': 'Zackular 2014, CRC',
         'crc_zeller': 'Zeller 2014, CRC',
         'edd_singh': 'Singh 2015, EDD',
         'hiv_dinh': 'Dinh 2015, HIV',
         'hiv_lozupone': 'Lozupone 2013, HIV',
         'hiv_noguerajulian': 'Noguera-Julian 2016, HIV',
         'ibd_gevers': 'Gevers 2014, IBD',
         'ibd_morgan': 'Morgan 2012, IBD',
         'ibd_papa': 'Papa 2012, IBD',
         'ibd_willing': 'Willing 2009, IBD',
         'mhe_zhang': 'Zhang 2013, LIV',
         'nash_wong': 'Wong 2013, NASH',
         'nash_zhu': 'Zhu 2013, NASH',
         'ob_goodrich': 'Goodrich 2014, OB',
         'ob_ross': 'Ross 2015, OB',
         'ob_turnbaugh': 'Turnbaugh 2009, OB',
         'ob_zhu': 'Zhu 2013, OB',
         'ob_zupancic': 'Zupancic 2012, OB',
         'par_scheperjans': 'Scheperjans 2015, PAR',
         'ra_scher': 'Scher 2013, ART',
         't1d_alkanani': 'Alkanani 2015, T1D',
         't1d_mejialeon': 'Mejia-Leon 2014, T1D'}
    return {i: d[i] for i in dataset_order}

def get_labeldict_for_overlap(dataset_order):
    """Dictionary for dataset labels for the percent overlap figure."""
    d = {'asd_kang': 'Kang (ASD)',
         'asd_son': 'Son (ASD)',
         'cdi_schubert': 'Schubert (CDI)',
         'cdi_singh': 'Singh (EDD)',
         'cdi_vincent': 'Vincent (CDI)',
         'cdi_youngster': 'Youngster (CDI)',
         'crc_baxter': 'Baxter (CRC)',
         'crc_chen': 'Chen (CRC)',
         'crc_wang': 'Wang (CRC)',
         'crc_zackular': 'Zackular (CRC)',
         'crc_zeller': 'Zeller (CRC)',
         'edd_singh': 'Singh (EDD)',
         'hiv_dinh': 'Dinh (HIV)',
         'hiv_lozupone': 'Lozupone (HIV)',
         'hiv_noguerajulian': 'Noguera-Julian (HIV)',
         'ibd_gevers': 'Gevers (IBD)',
         'ibd_morgan': 'Morgan (IBD)',
         'ibd_papa': 'Papa (IBD)',
         'ibd_willing': 'Willing (IBD)',
         'mhe_zhang': 'Zhang (LIV)',
         'nash_wong': 'Wong (NASH)',
         'nash_zhu': 'Zhu (NASH)',
         'ob_goodrich': 'Goodrich (OB)',
         'ob_ross': 'Ross (OB)',
         'ob_turnbaugh': 'Turnbaugh (OB)',
         'ob_zhu': 'Zhu (OB)',
         'ob_zupancic': 'Zupancic (OB)',
         'par_scheperjans': 'Scheperjans (PAR)',
         'ra_scher': 'Scher (ART)',
         't1d_alkanani': 'Alkanani (T1D)',
         't1d_mejialeon': 'Mejia-Leon (T1D)'}
    return {i: d[i] for i in dataset_order}

def get_phylo_colors(keep_rows):
    """
    Return the df with phylum and orders. Make color_dict to plot bar phylogeny.

    Parameters
    ----------
    keep_rows : list or pandas Index
        list of OTUs to return colors for. Should be the full taxxonomy
        (starting with k__;...) and in the order that you'll want to plot them.

    Return
    ------
    phylodf : pandas DataFrame
        dataframe with 'phylum', 'class', 'order', etc columns, in the same
        order is in keep_rows. Values are text like ('o__Clostridiales')
    color_dict : dict
        dictionary with {tax_level: RGBA tuple}, where tax_level
        is a string like 'o__Clostridiales'. All phyla and orders
        in keep_rows are in this dict.
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
    # orange (fuso), yellow, brown (eury), pink (verruco),
    # gray (teneri, cyano, lenti, synerg)
    # this color dict is for 'Set1' color palette
    # This first definition was for use with a ListedColorMap...
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
    # Convert directly to RGB values
    color_dict = {i: colors[color_dict[i]] for i in color_dict}

    ## Check that all phyla in phylodf have a color
    missingphyla = [i for i in phylodf['phylum'].unique()
                    if i not in color_dict]
    if len(missingphyla) > 0:
        print('You need to give the following phyla colors in get_phylo_colors():')
        print('\n'.join(missingphyla))
        raise ValueError

    order_color_dict = {}
    for p, subdf in phylodf.groupby('phylum'):
        # Get unique orders in that phylum
        orders = subdf['order'].unique()
        # Seed a light palette from the top-level phylum color.
        order_colors = sns.light_palette(color_dict[p], len(orders) + 2)[1:-1]
        # Assign corresponding index in `colors` to each order,
        # and add to the master color_dict
        color_dict.update(
            {orders[i]: order_colors[i] for i in range(len(orders))})

    return phylodf, color_dict

### Table writing stuff
def convert_to_latex(row):
    return ' & '.join(['\%'.join('\_'.join(i.split('_')).split('%')) for i in row.astype(str)]) + ' \\\ '

def write_latex_table(df, outfile):
    """
    Write the table in df as a latex table. Don't include any column
    headers (assuming that the user will do these manually in their
    manuscript...)
    """
    latex_series = df.apply(convert_to_latex, axis=1)
    with open(outfile, 'w') as f:
        f.write('\n'.join(latex_series))
    return None

def write_markdown_table(df, outfile):
    """
    Write the table in df as a markdown table. *Do* include the column
    names because github is less formal than manuscripts. :)
    """
    github_str = df.apply(lambda row: ' | '.join([i for i in row.astype(str)]), axis=1)
    with open(outfile, 'w') as f:
        f.write(' | '.join(df.columns) + '\n')
        f.write('|'.join(['----------']*len(df.columns)) + '\n')
        f.write('\n'.join(github_str) + '\n')
    return None
