#!/usr/bin/env python
"""
This script reads in a file with all qvalues for all datasets
and writes out a file with the meta-analysis results.
"""
import os
import argparse
import pandas as pd
import numpy as np

def count_sig(allresults, qthresh=0.05):
    """
    Count how often bacteria are significant in each disease.

    Parameters
    ---------
    allresults : pandas dataframe
        datasets in columns, genera in rows, signed q-values in matrix (signed according to effect direction)
    qthresh : float
        significance threshold

    Returns
    -------
    meta_counts : pandas DataFrame
        Dataframe indicating how often each genus is significant in each
        diseases. Has columns ['otu', 'disease', 'significant',
        'num_times_sig', 'genus'], where:
            'significant' is [-1, 1]
            'num_time_sigs' is the number of times each otu is significant
                in each disease/direction combination
            'otu' is the full OTU name, 'genus' is just the genus
    """
    def thresh_map(x):
        # Note: in upstream steps, pvalues of 0 were converted to 1e-20
        # if x is zero because the effect was zero, np.sign(x) returns 0. So we're good.
        if abs(x) <= qthresh:
            return np.sign(x)
        elif abs(x) > qthresh:
            return 0
        else:
            return x

    meta_results = allresults.applymap(thresh_map)
    meta_results['otu'] = meta_results.index
    meta_results = pd.melt(meta_results, id_vars='otu',
                           var_name='dataset', value_name='significant')

    # Replace edd_singh with cdi_singh
    meta_results = meta_results\
        .replace('edd_singh', 'cdi_singh')\
        .replace('noncdi_schubert', 'cdi_schubert2')
    meta_results['disease'] = meta_results['dataset']\
        .apply(lambda x: x.split('_')[0])
    # Drop rows with either nan or 0 in 'significant' column (i.e. not significant, no effect)
    meta_results = meta_results.dropna()
    meta_results = meta_results.loc[meta_results['significant'].isin([-1,1])]

    # Get number of times each OTU is significant in each disease
    # (for each direction)
    meta_counts = meta_results.groupby(['otu', 'disease', 'significant']).size()
    meta_counts.name = 'num_times_sig'
    meta_counts = meta_counts.reset_index()
    meta_counts['genus'] = meta_counts['otu'].apply(lambda x: x.split(';')[-1])

    return meta_counts

def within_disease_meta_analysis(meta_counts, all_otus=None,
                                 diseases=['cdi', 'ibd', 'crc', 'ob', 'hiv']):
    """
    Find meta-significant bugs for given diseases, based on counting how often
    they're significant within specified diseases.

    Parameters
    ----------
    meta_counts : pandas DataFrame
        dataframe with columns 'disease', 'significant', and 'num_times_sig'
    all_otus : list or pandas Index
        list of otus to use as index of disease_df. Note: meta_counts
        doesn't have all the original OTUs (from the qvalues dataframe),
        so if you want a result with the same shape as the original qvalues,
        you need to specify the OTUs/rows here. Default is all unique OTUs
        in meta_counts (i.e. just the ones that were significant in at least
        one study)
    diseases : list
        list of diseases to perform meta-analysis for

    Returns
    -------
    disease_df : pandas dataframe
        Dataframe with diseases in columns and genera in rows.
        Values indicate whether that genus is disease-significant.
        Genera are disease-significant if they are significant
        in at least 2 studies within each disease (in the same
        direction).
    """

    ## For each of my diseases of interest, get the genera
    ## which are significant in at least 2 studies
    meta_bugs = {}
    for dis in diseases:
        meta_bugs[dis] = dict.fromkeys(['healthy', 'disease'])
        meta_bugs[dis]['healthy'] = meta_counts \
            .query('disease == @dis') \
            .query('significant == -1') \
            .query('num_times_sig >= 2')['otu'].tolist()
        meta_bugs[dis]['disease'] = meta_counts \
            .query('disease == @dis') \
            .query('significant == 1') \
            .query('num_times_sig >= 2')['otu'].tolist()

    # Put these 'disease-specific' significant bugs into dataframe with -1/1 values
    if all_otus is None:
        all_otus = meta_counts['otu'].unique()

    disease_df = pd.DataFrame(index=all_otus,
                              columns=meta_bugs.keys())
    for dis in meta_bugs:
        purple_bugs = [i for i in meta_bugs[dis]['healthy']
            if i in meta_bugs[dis]['disease']]
        if len(purple_bugs) > 1:
            print('Warning!! Some disease-specific bugs are significant in both directions!! Go check this!!')
        disease_df.loc[meta_bugs[dis]['healthy'], dis] = -1
        disease_df.loc[meta_bugs[dis]['disease'], dis] = 1
    disease_df = disease_df.astype(float)

    return disease_df

def cross_disease_meta_analysis(meta_counts, num_diseases=2, exclude_dis=None,
                                all_otus=None):
    """
    Get core bacteria which are significant in the same direction in at least
    `num_diseases` different diseases.

    Parameters
    ----------
    meta_counts : pandas DataFrame
        dataframe with columns 'disease', 'significant', and 'num_times_sig'
    num_diseases : int
        number of diseases a genus must be significant in
        to be considered a "core" bug. Default is 2 (i.e. must
        be significant at least once in at least 2 diseases)
    exclude_dis : str or list of strings
        diseases to exclude from the calculation, as in 'disease' column of
        meta_counts.
    all_otus : list or pandas Index
        list of otus to use as index of disease_df. Note: meta_counts
        doesn't have all the original OTUs (from the qvalues dataframe),
        so if you want a result with the same shape as the original qvalues,
        you need to specify the OTUs/rows here. Default is all unique OTUs
        in meta_counts (i.e. just the ones that were significant in at least
        one study)

    Returns
    -------
    overall_df : pandas dataframe
        Dataframe has one column, 'overall'.
        Values indicate overall significance: 1 (disease),
        -1 (health), and 0 (both).
        Genera are considered overall significant if they are
        significant in the same direction in at least
        n_diseases *diseases* (regardless of number of studies
        within each disease).
        If a genus is overall significant in both health- and
        disease-directions, then it's given a value of 0.
        In other words, if it's associated with disease in
        at least n_diseases diseases *and also* associated with
        health in at least n_diseases diseases, it gets a 0.
    """

    if exclude_dis is not None:
        meta_counts = meta_counts.query('disease != @exclude_dis')

    ## Get cross-disease bugs:
    # Find bugs which are significant in at least n_diseases
    cross_dis_counts = []
    # Get number of times each OTU is significant in each dataset (for each direction)
    # meta_counts is already grouped by otu, disease, and significance direction
    for o, subdf in meta_counts.groupby('otu'):
        num_h = subdf.query('significant == -1').shape[0]
        num_d = subdf.query('significant == 1').shape[0]
        cross_dis_counts.append([o, num_h, num_d])
    cross_dis_counts = pd.DataFrame(data=cross_dis_counts,
        columns=['otu', 'dis_healthy', 'dis_disease'])

    overall_mixed = cross_dis_counts \
        .query('dis_healthy >= @num_diseases') \
        .query('dis_disease >= @num_diseases') \
        ['otu'].tolist()
    overall_health = [i for i in
        cross_dis_counts.query('dis_healthy >= @num_diseases')['otu'].tolist()
        if i not in overall_mixed]
    overall_disease = [i for i in
        cross_dis_counts.query('dis_disease >= @num_diseases')['otu'].tolist()
        if i not in overall_mixed]

    # Put "overall" significant bugs into dataframe with -1/1/0 values
    if all_otus is None:
        all_otus = meta_counts['otu'].unique()
    overall_df = pd.DataFrame(index=all_otus, columns=['overall'])
    overall_df.loc[overall_health] = -1
    overall_df.loc[overall_disease] = 1
    overall_df.loc[overall_mixed] = 0
    overall_df = overall_df.astype(float)

    return overall_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('qvalues', help='file with qvalues; genera in rows, datasets in columns')
    parser.add_argument('out_dir', help='directory to save results files to')
    parser.add_argument('qthresh', help='significance threshold', default=0.05, type=float)
    parser.add_argument('n_diseases', help='number of diseases to use in calculating "core" genera', default=2, type=int)
    parser.add_argument('--no-cdi', help='flag to exclude diarrhea datasets '
        + 'in determining core bugs.', action='store_true')
    parser.add_argument('--exclude-nonhealthy', help='flag to exclude '
        + 'studies without healthy controls and hiv_lozupone from the '
        + 'overall cross-disease meta-analysis', action='store_true')
    parser.add_argument('--disease', help='flag to perform disease-wise '
        + 'meta-analysis.', action='store_true')
    parser.add_argument('--overall', help='flag to perform cross-disease '
        + 'meta-analysis.', action='store_true')

    args = parser.parse_args()

    qvals = pd.read_csv(args.qvalues, sep='\t', index_col=0)

    meta_counts = count_sig(qvals, args.qthresh)

    # Disease-specific bugs
    if args.disease:
        disease_df = within_disease_meta_analysis(
            meta_counts, all_otus=qvals.index)

        # Save file
        disease_out = os.path.join(args.out_dir,
            'meta.counting.q-{}.disease_wise.txt'.format(args.qthresh))
        disease_df.to_csv(disease_out, sep='\t')

    if args.overall:
        ## Re-count how many times each bug is significant in each disease, after
        ## excluding any datasets without healthy controls.
        if args.exclude_nonhealthy:
            to_exclude = ['ibd_papa', 'ibd_gevers', 'hiv_lozupone']
            qvals = qvals.drop(to_exclude, axis=1)
            meta_counts = count_sig(qvals, args.qthresh)

        # Core bugs
        if args.no_cdi:
            overall_df = cross_disease_meta_analysis(
                meta_counts, args.n_diseases, exclude_dis=['cdi'])
        else:
            overall_df = cross_disease_meta_analysis(
                meta_counts, args.n_diseases)

        if args.no_cdi:
            overall_out = os.path.join(args.out_dir,
                'meta.counting.q-{}.{}_diseases.across_all_diseases_except_cdi.txt'\
                    .format(
                        args.qthresh, args.n_diseases))
        else:
            overall_out = os.path.join(args.out_dir,
                'meta.counting.q-{}.{}_diseases.across_all_diseases.txt'.format(
                    args.qthresh, args.n_diseases))

        overall_df.to_csv(overall_out, sep='\t')
