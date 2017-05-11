#!/usr/bin/env python
"""
This script calculates the ubiquity and abundance of genera across
healthy, case, and all patients. It writes a tidy dataframe.
"""
import argparse
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('TKAgg')
import seaborn as sns

# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
import FileIO as fio
import util

def read_all_and_return_abun_ubiquity(datadir, fnpvals):
    """
    Read all clean datasets in datadir and return a tidy dataframe
    with the various ubiquity/abundance calculations for each genus.

    Parameters
    ----------
    datadir : str
        path to directory with *.otu_table.clean.feather and
        *.metadata.clean.feather files
    fnpvals : str
        path to file with 'overall' significant bugs (should have column labeled
        'overall' and OTUs in rows)

    Returns
    -------
    tidy : pandas tidy dataframe
       has the following columns:
          otu: str, 'k__Bacteria;...g__Akkermansia'
          variable: str, 'abundance_from_pooled_mean_total',
                    'ubiquity_mean_of_datasets_h', etc...
          value: float
          metric: str, {'abundance', 'ubiquity'}
          calculation: str, {'from_pooled_mean', 'mean_of_datasets'}
          patient: {'total', 'dis', 'h'}
          overall: float or str {'nan', 1.0, -1.0, 0.0}
          color: float or str, {RGBA values or 'k'}
          overall_significance: str, {'not_sig', 'disease', 'health', 'mixed'}
    """
    toconcat = []
    print('Reading datasets...')
    datasetids = fio.get_dataset_ids(datadir)
    for dataset in datasetids:
        print(dataset),
        ## Read dataset
        df, meta = fio.read_dataset_files(dataset, datadir)
        df = util.raw2abun(df)

        ## Collapse to genus level
        df = util.collapse_taxonomic_contents_df(df, 'genus')
        classes_list = fio.get_classes(meta)
        [ctrl_smpls, dis_smpls] = fio.get_samples(meta, classes_list)

        ## Do calculations for all patients
        # Note: dis_smpls and ctrl_smpls sometimes just grabs a subset of
        # patients. e.g. in CRC studies, this discards adenoma patients
        #all_smpls = ctrl_smpls + dis_smpls
        all_smpls = list(df.index)
        # tmp is a DataFrame with the calculations for different patient
        # groups in the columns, genera in rows
        tmp = pd.DataFrame(df.loc[all_smpls].sum(),
                           columns=['total_total_abun'])
        tmp['total_present'] = df.loc[all_smpls]\
                                 .applymap(lambda x: 1 if x else 0)\
                                 .sum()
        tmp['total_samples'] = len(all_smpls)

        ## Do calculations across healthy patients only
        H_smpls = meta.query('DiseaseState == "H"').index
        if len(H_smpls) > 0:
            tmp['total_h_abun'] = df.loc[H_smpls].sum()
            tmp['h_present'] = df.loc[H_smpls]\
                                 .applymap(lambda x: 1 if x else 0)\
                                 .sum()
            tmp['h_samples'] = len(H_smpls)

        else:
            # add the non-healthy controls to our disease patients
            dis_smpls += ctrl_smpls

        ## Do calculations across non-healthy patients
        tmp['total_dis_abun'] = df.loc[dis_smpls].sum()
        tmp['dis_present'] = df.loc[dis_smpls]\
                               .applymap(lambda x: 1 if x else 0)\
                               .sum()
        tmp['dis_samples'] = len(dis_smpls)

        tmp['dataset'] = dataset
        tmp['otu'] = tmp.index

        toconcat.append(tmp)

    df = pd.concat(toconcat, ignore_index=True)

    # Calculate ubiquity and abundance metrics for genera
    # df now has many columns, see the calculate_ubiquity_and_abun docstring
    # for these columns
    df = calculate_ubiquity_and_abun(df, 'h')
    df = calculate_ubiquity_and_abun(df, 'dis')
    df = calculate_ubiquity_and_abun(df, 'total')

    # Turn into tidy df with 'metric', 'calculation', and 'patient' columns
    tidy = tidyfy_df(df)

    # Add the core significance value for each genus
    overallsig = pd.read_csv(fnpvals, sep='\t', index_col=0)

    tidy = tidy.merge(overallsig, right_index=True, left_on='otu', how='left')
    # replace NaN with 'nan' to allow dict hashing
    tidy['overall'] = tidy['overall'].replace(np.nan, 'nan')
    # color each genus by its overall status (not used in final plots)
    pal = sns.diverging_palette(220,20, n=7)
    colordict = {1: tuple(pal[-1]),
                -1: tuple(pal[0]),
                0: 'k',
                'nan': tuple(pal[3])}
    tidy['color'] = tidy['overall'].map(lambda x: colordict[x])
    # add human-readable 'overall' labels
    rename_overall = {1: 'disease', -1: 'health', 'nan': 'not_sig', 0: 'mixed'}
    tidy.loc[tidy.index, 'overall_significance'] = \
        tidy['overall'].map(lambda x: rename_overall[x])

    return tidy

def calculate_ubiquity_and_abun(df, patient_type):
    """
    Use groupby-fu to calculate the ubiquity and abundance of genera
    across patients of patient_type.

    Parameters
    ----------
    df : pandas dataframe
        df with the following columns: '<H/dis/total>_present',
        '<H/dis/total>_samples', 'total_<H/dis/total>_abun', and 'otu'
    patient_type : str
        'H', 'dis', or 'total'. Type of patient to calculate stuff over.

    Returns
    -------
    df : pandas dataframe
        modified input df (not modified in place) with additional columns:
        'ubiquity_in_one_dataset_<h/dis/total>',
        'mean_abundance_in_one_dataset_<h/dis/total>',
        'ubiquity_mean_of_datasets_<h/dis/total>',
        'ubiquity_from_pooled_calc_<h/dis/total>',
        'abundance_mean_of_dataset_<h/dis/total>',
        'abundance_from_pooled_calc_<h/dis/total>'
    """
    if patient_type == 'H' or patient_type == 'h':
        coltype = 'h'
    elif patient_type == 'dis':
        coltype = 'dis'
    elif patient_type == 'total':
        coltype = 'total'
    else:
        print('Unknown patient type.')

    ## Calculate various metrics per OTU using different types of averages
    #(average of dataset-wise average, overall average)

    ## Ubiquity of OTUs
    # samples with OTU present / total samples
    df['ubiquity_in_one_dataset_' + coltype] = \
        df[coltype + '_present'].astype(float)/df[coltype + '_samples']

    # Get the mean of individual dataset's ubiquities, for each OTU
    # Turn that series into a dataframe and merge it with the existing df
    # on the 'otu' column to add the 'ubiquity_mean_of_datasets' column
    # to the original df
    df = (df.groupby('otu') \
            .mean()['ubiquity_in_one_dataset_' + coltype] \
            .to_frame('ubiquity_mean_of_datasets_' + coltype)) \
            .merge(df,right_on='otu', left_index=True)

    # Overall ubiquity is the sum of all samples with the OTU present,
    # divided by all the samples total
    df = (df.groupby('otu').sum()[coltype + '_present'] \
          / df.groupby('otu').sum()[coltype + '_samples']) \
          .to_frame(name='ubiquity_from_pooled_mean_' + coltype) \
          .merge(df, right_on='otu', left_index=True)

    ## Abundance - consider abundance across *only* samples who have the bug
    # Mean of dataset-wise mean abundances
    # First get the mean abundance of the OTU in each dataset
    df['mean_abundance_in_one_dataset_' + coltype] = \
            df['total_' + coltype + '_abun'] / df[coltype + '_present']

    # Then take the mean of those means for each OTU
    # And merge with the original df to create 'abundance_mean_of_datasets' col
    df = (df.groupby('otu') \
            .mean()['mean_abundance_in_one_dataset_' + coltype] \
            .to_frame('abundance_mean_of_datasets_' + coltype)) \
            .merge(df,right_on='otu', left_index=True)

    # Mean abundance overall: total abundance across all studies divided
    # by total number of samples with the OTU present across all studies
    df = (df.groupby('otu').sum()['total_' + coltype + '_abun'] \
          / df.groupby('otu').sum()[coltype + '_present']) \
          .to_frame(name='abundance_from_pooled_mean_' + coltype) \
          .merge(df, right_on='otu', left_index=True)

    return df

def tidyfy_df(df):
    """
    Returns tidy df pivoted around 'otu', for the aggregate ubiquity and abundance
    measures.

    Input df should have columns labeled like "ubiquity_calc_type_patients" or
    "abundance_calc_type_patients" where the first underscore-delimited value
    is "abundance" or "ubiquity" and the last one is "dis", "h", or "total"
    (or some other patient type indicator). The middle values are the type of
    calculation used (e.g. "from_pooled_calc", "mean_of_datasets")

    Note that columns with 'in_one_dataset' are discarded.
    """

    id_vars = ['otu']
    value_vars = [i for i in df.columns if i.startswith('ubiquity') or i.startswith('abundance')]
    value_vars = [i for i in value_vars if 'in_one_dataset' not in i]

    tidydf = pd.melt(df, id_vars=id_vars, value_vars=value_vars).drop_duplicates()

    tidydf['metric'] = tidydf['variable'].apply(lambda x: x.split('_')[0])
    tidydf['calculation'] = tidydf['variable'].apply(lambda x: x.split('_',1)[1].rsplit('_',1)[0])
    tidydf['patient'] = tidydf['variable'].apply(lambda x: x.split('_')[-1])

    return tidydf

p = argparse.ArgumentParser(description="This script calculates the ubiquity "
    "and abundance of genera across healthy, case, and all patients. It "
    "writes a tidy dataframe.")
p.add_argument('datadir', help='directory with clean OTU and metadata tables')
p.add_argument('fnoverall', help='file with the "overall" significance')
p.add_argument('out', help='file to write results to')

args = p.parse_args()

tidy = read_all_and_return_abun_ubiquity(args.datadir, args.fnoverall)
tidy.to_csv(args.out, sep='\t')
