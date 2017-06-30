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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('clean_data_dir', help='directory with clean OTU and metadata tables')
    parser.add_argument('out_dir', help='directory to save output q values file to')
    args = parser.parse_args()

    qthresh = 0.05
    stats_method = 'kruskal-wallis'

    dfdict = fio.read_dfdict_data(args.clean_data_dir)

    print('Doing univariate tests...')
    for dataset in dfdict:
        df = dfdict[dataset]['df']
        meta = dfdict[dataset]['meta']

        # Collapse to genus level
        df = util.collapse_taxonomic_contents_df(df, 'genus')

        # Get samples in each class
        classes_list = fio.get_classes(meta)
        if len(classes_list[0]) == 0 or len(classes_list[1]) == 0:
            raise ValueError('Something wrong with ' + dataset + ' metadata.')

        # Initialize results sub-dict
        dfdict[dataset]['comparisons'] = {}

        for dis_label in classes_list[1]:
            comparison = classes_list[0][0] + '_vs_' + dis_label

            sub_list = [classes_list[0], [dis_label]]
            H_smpls, dis_smpls = fio.get_samples(meta, sub_list)

            # Do some stats. 'results' is a dataframe with two
            # columns, 'p' and 'test-stat'
            results = util.compare_otus_teststat(
                df, H_smpls, dis_smpls, method=stats_method,
                multi_comp='fdr')

            dfdict[dataset]['comparisons'][comparison] = \
                {'df': df, 'dis_smpls': dis_smpls,
                 'H_smpls': H_smpls, 'results': results}

    ## Manipulate each dataset's results with values signed
    ## according to median and mean effects in dis - H
    med_allresults_lst = None
    mean_allresults_lst = None
    allpvals = None
    for dataset in dfdict:
        for comp in dfdict[dataset]['comparisons']:
            ## Manipulate the results dataframe
            # into a df with signed values according to direction of change
            # Calculates using both mean and median.
            results = dfdict[dataset]['comparisons'][comp]['results']
            df = dfdict[dataset]['comparisons'][comp]['df']
            dis_smpls = dfdict[dataset]['comparisons'][comp]['dis_smpls']
            H_smpls = dfdict[dataset]['comparisons'][comp]['H_smpls']

            med_results, mean_results = sign_results(results, df, dis_smpls, H_smpls, dataset, col='q')

            # Update df names to have both the dataset and the comparison
            med_results.name = dataset + '-' + comp
            mean_results.name = dataset + '-' + comp

            try:
                med_allresults_lst.append(med_results)
                mean_allresults_lst.append(mean_results)
            except AttributeError:
                med_allresults_lst = [med_results]
                mean_allresults_lst = [mean_results]

    print('Doing univariate tests... Finished')

    ## Concat list of signed results series into a dataframe with genera in rows, datasets in columns
    med_allresults =  pd.concat(med_allresults_lst, axis=1)
    mean_allresults =  pd.concat(mean_allresults_lst, axis=1)

    medfile = os.path.join(args.out_dir, 'q-val_all_results.median.{}.case-control.txt'.format(stats_method))
    meanfile = os.path.join(args.out_dir, 'q-val_all_results.mean.{}.case-control.txt'.format(stats_method))
    med_allresults.to_csv(medfile, sep='\t')
    mean_allresults.to_csv(meanfile, sep='\t')
