#!/usr/bin/env python
"""
This script calculates q-values for all genera in all datasets.
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd

# Add this repo to the path
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
import util
import FileIO as fio


def reformat_results(results, col):
    """
    Takes dataframe results with columns 'col' and 'effect'
    Returns +/- values in 'col', signed according to 'effect'
    """

    ## Need to convert any zeros to a very small number, so that
    ## number*np.sign(effect) is non-zero.
    # Sanity check: if my pvalue is 0, and I replace that with
    # 1e-20, none of my downstream analyses will be affected.
    # Also, downstream code (clean_results_df) removes rows with all zero's, i.e. where the effect was zero.
    # So this is double-good, bc now pvalues = 0 won't be removed later.
    results[col] = results[col].replace(0, 1e-20)
    return results[col]*np.sign(results['effect'])

def sign_results(results, df, dis_smpls, H_smpls, dataset, col):
    """
    Convert one results dataframe with p-values for all OTUs in df
    into two dataframes with 'effect' column calculated using median and mean.

    Parameters
    ----------
    results : pandas dataframe
        Dataframe with genera in index and column 'col'
    df : pandas dataframe
        Dataframe with genera in columns and samples in index
    dis, H_smpls : lists
        Lists of samples in each class.
        Samples should be in the index of df.
    dataset : str
        Dataset ID. Used to name the returned Series.
    col : str
        Which column in results dataframe to return signed values for

    Returns
    -------
    med/mean_results : pandas Series
        Series with +/- 'col' for each genus.
        Sign corresponds to direction of effect, calculated both by mean and median.
        Positive indicates higher in disease, negatives is higher in healthy.
    """
    # Median results
    results['effect'] = df.loc[dis_smpls].median() - df.loc[H_smpls].median()
    med_results = reformat_results(results, col)
    med_results.name = dataset

    # Mean results
    results['effect'] = df.loc[dis_smpls].mean() - df.loc[H_smpls].mean()
    mean_results = reformat_results(results, col)
    mean_results.name = dataset

    return med_results, mean_results

parser = argparse.ArgumentParser()
parser.add_argument('clean_data_dir', help='directory with clean OTU and '
    + ' metadata tables')
parser.add_argument('out_file', help='path to write qvalues to')
parser.add_argument('--subset', help='file with list of dataset IDs that '
    + 'should be analyzed, if you want to analyze only a subset of the '
    + 'datasets in clean_data_dir. Dataset IDs should match IDs in file '
    + 'names, and be one per line.', default=None)
parser.add_argument('--split-cases', help='flag to analyze each case type '
    + 'separately.', action='store_true')
args = parser.parse_args()

qthresh = 0.05
stats_method = 'kruskal-wallis'

dfdict = fio.read_dfdict_data(args.clean_data_dir, subset=args.subset)

print('Doing univariate tests...')
resultsdict = {}

for dataset in dfdict:
    df = dfdict[dataset]['df']
    meta = dfdict[dataset]['meta']

    # Collapse to genus level
    df = util.collapse_taxonomic_contents_df(df, 'genus')


    if args.split_cases:
        # Get samples in each class. Note that the two fio.get_classes
        # functions are basically the same, just with different diseases
        # hard-coded in.
        classes_list = fio.get_classes(meta, dataset)

        # Go through each case group one by one
        for dis_label in classes_list[1]:
            # old dataset = ibd_alm, new dataset = uc_alm
            newdataset = dis_label.lower() + '_' + dataset.split('_')[1]

            # Get samples
            sub_list = [classes_list[0], [dis_label]]
            H_smpls, dis_smpls = fio.get_samples(meta, sub_list)

            # Do some stats. 'results' is a dataframe with two
            # columns, 'p' and 'test-stat'
            results = util.compare_otus_teststat(
                df, H_smpls, dis_smpls, method=stats_method,
                multi_comp='fdr')

            resultsdict[newdataset] = {
                 'df': df, 'dis_smpls': dis_smpls,
                 'H_smpls': H_smpls, 'results': results}

    else:
        classes_list = fio.get_classes(meta, dataset)
        # Just combine all cases together
        H_smpls, dis_smpls = fio.get_samples(meta, classes_list)

        # Do some stats. 'results' is a dataframe with two
        # columns, 'p' and 'test-stat'
        results = util.compare_otus_teststat(
            df, H_smpls, dis_smpls, method=stats_method,
            multi_comp='fdr')

        resultsdict[dataset] = {
            'df': df, 'meta': meta,
            'dis_smpls': dis_smpls, 'H_smpls': H_smpls,
            'results': results}

## Manipulate each dataset's results with values signed
## according to median and mean effects in dis - H
med_allresults_lst = []
mean_allresults_lst = []

for dataset in resultsdict:
    ## Manipulate the results dataframe
    # into a df with signed values according to direction of change
    # Calculates using both mean and median.
    results = resultsdict[dataset]['results']
    df = resultsdict[dataset]['df']
    dis_smpls = resultsdict[dataset]['dis_smpls']
    H_smpls = resultsdict[dataset]['H_smpls']

    med_results, mean_results = \
        sign_results(results, df, dis_smpls, H_smpls, dataset, col='q')

    med_allresults_lst.append(med_results)
    mean_allresults_lst.append(mean_results)

print('Doing univariate tests... Finished')

## Concat list of signed results series into a dataframe with
## genera in rows, datasets in columns
med_allresults =  pd.concat(med_allresults_lst, axis=1)
mean_allresults =  pd.concat(mean_allresults_lst, axis=1)

# Note: I don't ever use the median-based results in the paper, but this
# is where you could find them.
#med_allresults.to_csv(medfile, sep='\t')
mean_allresults.to_csv(args.out_file, sep='\t')
