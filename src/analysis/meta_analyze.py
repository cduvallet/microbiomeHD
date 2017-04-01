#!/usr/bin/env python
"""
This script reads in a file with all qvalues for all datasets
and writes out a file with the meta-analysis results.
"""
import os
import argparse
import pandas as pd
import numpy as np

def count_based_meta_analysis(allresults, qthresh,
                              diseases=['cdi', 'ibd', 'crc', 'ob', 'hiv'],
                              num_diseases=2):
    """
    Find meta-significant bugs based on counting how often they're significant within and across diseases.

    Parameters
    ---------
    allresults : pandas dataframe
        datasets in columns, genera in rows, signed q-values in matrix (signed according to effect direction)
    qthresh : float
        significance threshold
    diseases : list
        list of diseases to perform meta-analysis on
    num_diseases : int
        number of diseases a genus must be significant in
        to be considered a "core" bug. Default is 2 (i.e. must
        be significant at least once in at least 2 diseases)

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

    disease_df : pandas dataframe
        Dataframe with diseases in columns and genera in rows.
        Values indicate whether that genus is disease-significant.
        Genera are disease-significant if they are significant
        in at least 2 studies within each disease (in the same
        direction).
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
    meta_results = pd.melt(meta_results, id_vars='otu', var_name='dataset', value_name='significant')

    # Replace edd_singh with cdi_singh
    meta_results = meta_results.replace('edd_singh', 'cdi_singh')
    meta_results['disease'] = meta_results['dataset'].apply(lambda x: x.split('_')[0])
    # Drop rows with either nan or 0 in 'significant' column (i.e. not significant, no effect)
    meta_results = meta_results.dropna()
    meta_results = meta_results.loc[meta_results['significant'].isin([-1,1])]

    # Get number of times each OTU is significant in each disease (for each direction)
    meta_counts = meta_results.groupby(['otu', 'disease', 'significant']).size()
    meta_counts.name = 'num_times_sig'
    meta_counts = meta_counts.reset_index()
    meta_counts['genus'] = meta_counts['otu'].apply(lambda x: x.split(';')[-1])

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
    disease_df = pd.DataFrame(index=allresults.index,
                              columns=meta_bugs.keys())
    for dis in meta_bugs:
        purple_bugs = [i for i in meta_bugs[dis]['healthy'] if i in meta_bugs[dis]['disease']]
        if len(purple_bugs) > 1:
            print('Warning!! Some disease-specific bugs are significant in both directions!! Go check this!!')
        disease_df.loc[meta_bugs[dis]['healthy'], dis] = -1
        disease_df.loc[meta_bugs[dis]['disease'], dis] = 1
    disease_df = disease_df.astype(float)

    ## Get cross-disease bugs:
    # Find bugs which are significant in at least n_diseases
    cross_dis_counts = []
    # Get number of times each OTU is significant in each dataset (for each direction)
    # meta_counts is already grouped by otu, disease, and significance direction
    for o, subdf in meta_counts.groupby('otu'):
        num_h = subdf.query('significant == -1').shape[0]
        num_d = subdf.query('significant == 1').shape[0]
        cross_dis_counts.append([o, num_h, num_d])
    cross_dis_counts = pd.DataFrame(data=cross_dis_counts, columns=['otu', 'dis_healthy', 'dis_disease'])

    overall_mixed = cross_dis_counts.query('(dis_healthy >= @num_diseases) & (dis_disease >= @num_diseases)')['otu'].tolist()
    overall_health = [i for i in cross_dis_counts.query('dis_healthy >= @num_diseases')['otu'].tolist() if i not in overall_mixed]
    overall_disease = [i for i in cross_dis_counts.query('dis_disease >= @num_diseases')['otu'].tolist() if i not in overall_mixed]

    # Put "overall" significant bugs into dataframe with -1/1/0 values
    overall_df = pd.DataFrame(index=allresults.index, columns=['overall'])
    overall_df.loc[overall_health] = -1
    overall_df.loc[overall_disease] = 1
    overall_df.loc[overall_mixed] = 0
    overall_df = overall_df.astype(float)

    return overall_df, disease_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('qvalues', help='file with qvalues; genera in rows, datasets in columns')
    parser.add_argument('out_dir', help='directory to save results files to')
    parser.add_argument('qthresh', help='significance threshold', default=0.05, type=float)
    parser.add_argument('n_diseases', help='number of diseases to use in calculating "core" genera', default=2, type=int)

    args = parser.parse_args()

    qvals = pd.read_csv(args.qvalues, sep='\t', index_col=0)

    overall_df, disease_df = count_based_meta_analysis(qvals, args.qthresh, num_diseases=args.n_diseases)

    overall_out = os.path.join(args.out_dir, 'meta.counting.q-{}.{}_diseases.across_all_diseases.txt'.format(args.qthresh, args.n_diseases))
    disease_out = os.path.join(args.out_dir, 'meta.counting.q-{}.disease_wise.txt'.format(args.qthresh))

    overall_df.to_csv(overall_out, sep='\t')
    disease_df.to_csv(disease_out, sep='\t')
