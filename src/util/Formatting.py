#!/usr/bin/env python
"""
Useful functions in formatting tables and figures.
"""
import numpy as np

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
            df['dataset'].apply(lambda x:  x.startswith(d)),
            'dataset'].values
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
