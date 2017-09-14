#!/usr/bin/env python
"""
This script looks at the expected concordance between all studies in
all pairwise disease comparisons.
"""
import argparse

import pandas as pd
import numpy as np

# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
from util import shuffle_col

def get_n_concord(df, col1, col2, intersection=False):
    """
    Get the number of concordant values in col1 and col2 of df.

    Values that are nan in one column and not nan in the other are
    considered discordant. Values that are nan in both columns are
    dropped.

    Note: if you've converted your dataframe of qvalues to +/- 1
    values, this code converts zeros to NaN's so that they are
    not considered concordant (np.nan == np.nan is False,
    but 0 == 0 is True).

    Parameters
    ----------
    df : pandas DataFrame
        df with genera in rows, studies in columns
    col1, col2 : str
    intersection : bool
        Whether to calculate concordance across all genera present in at least
        one of the columns (False) or only genera which are in *both*
        columns (True)
    """

    if intersection:
        tmp = df[[col1, col2]].dropna()
    else:
        tmp = df[[col1, col2]].dropna(how='all')

    tmp = tmp.replace(0, np.nan)

    # If neither column has anything significant
    if np.all(tmp.isnull()):
        return 0
    else:
        return sum(tmp[col1] == tmp[col2])

def get_pvals_df(df):
    """
    Calculate pvalues from the dataframe with the empirical distribution.

    Parameters
    ----------
    df : pandas DataFrame
        has columns ['study1', 'study2', 'qthresh', 'metric', 'n_concord',
        'dis1', 'dis2']. (Produced from concordance_analysis.py)
        The 'metric' column indicates whether the value in 'n_concord' is the
        actual value ("observed") or a shuffled, expected value ("expected").

    Returns
    -------
    pvalsdf : pandas DataFrame
        has columns ['study1', 'study2', 'p'] where 'p' is the p value of the
        expected n_concord being >= the observed n_concord.
        Note that pvalsdf contains *redundant* study1-study2 pairs (so that
        calling pvalsdf.pivot() yields a symmetric matrix)
    """

    # Make pvalues dataframe from the empirical distribution
    pvals = []
    for g, subdf in df.groupby(['study1','study2', 'qthresh']):
        n_obs = subdf.query('metric == "observed"')['n_concord'].values
        n_exp = subdf.query('metric == "expected"')['n_concord'].values
        p = sum(n_exp >= n_obs)/float(len(n_exp))

        pvals.append(
            list(g) + [p, subdf['dis1'].iloc[0], subdf['dis2'].iloc[0]])
    pvalsdf = pd.DataFrame(data=pvals,
        columns=['study1', 'study2', 'qthresh', 'p', 'dis1', 'dis2'])
    # Add study1-study2 flipped pairs so pivoted df is symmetric
    pvalsdf = pd.concat(
        (pvalsdf,
         pvalsdf[['study1', 'study2', 'p']].rename(
            columns={'study1': 'study2', 'study2': 'study1'})))\
        .drop_duplicates(['study1', 'study2', 'p'])
        # the within-disease comparisons get doubled here
    return pvalsdf

p = argparse.ArgumentParser()
p.add_argument('qvals', help='path to file with signed qvalues. Datasets '
    + 'in columns, genera in rows, signed qvalues in values.')
p.add_argument('--method', help='how to calculate pvalues. perm: '
    + 'permutation-based, p = sum(n_exp >= n_obs)/float(len(n_exp)). '
    + 'fisher_overlap: fisher exact test on genera which are present in '
    + 'both studies. fisher_all: fisher exact test on all genera, where '
    + 'missing genera are replaced with discordant directions. 3x3: '
    + 'chi squared test on distribution of {-1, 1, nan}, where genera '
    + 'which are missing in both studies are dropped [not implemented]. '
    + ' spearman: spearman correlation, automatically drops any NaNs. '
    + '[default: %(default)s]',
    choices=['perm', 'fisher', '3x3' 'cohen', 'spearman'],
    default='perm')
p.add_argument('--intersect', help='whether to calculate metric based on '
    + 'only genera which are present in both studies. Note that some '
    + 'metrics may not work without this... [default: %(default)s]',
    action='store_true')

p.add_argument('--qthresh', help='significance threshold for which genera to '
    + 'include in concordance analysis. qthresh = 1 means that effect '
    + 'directions will be compared (regardless of signifiance). qthresh = '
    + '0.05 means that directions of significant genera will be compared. '
    + '[default: %(default)s]', default=1.0, type=float)
p.add_argument('--nreps', help='number of shuffles to build empirical null. '
    + 'Only used if method == "perm". [default: %(default)s]',
    default=1000, type=int)
p.add_argument('--tidy_fout', help='file to write tidy results to. Only used '
    + 'if method is "perm".')
p.add_argument('pvals_fout', help='file to write pvalues to.')
args = p.parse_args()

# Read in qvalues
df = pd.read_csv(args.qvals, sep='\t', index_col=0)
# Rename edd_singh to cdi_singh
df = df.rename(columns={'edd_singh': 'cdi_singh',
                        'noncdi_schubert': 'cdi_schubert2'})


# Get diseases to compare between
diseases = list(set([i.split('_')[0] for i in df.columns]))

if args.method == 'perm':
    # Convert signed qvalues to +/- 1 if significant and 0 if not
    sigmap = lambda x: np.sign(x) if abs(x) <= args.qthresh else 0
    effdir = df.applymap(sigmap)

    dist = []
    ## Count number of actual concordant values
    print('Observed concordance...')
    for dis1 in diseases:
        dis1_studies = [i for i in df.columns if i.startswith(dis1)]
        for dis2 in diseases[diseases.index(dis1):]:
            dis2_studies = [i for i in df.columns if i.startswith(dis2)]

            for study1 in dis1_studies:
                for study2 in dis2_studies:
                    if study1 != study2:
                        n_concord = get_n_concord(effdir, study1, study2,
                            args.intersect)
                        dist.append(
                            [dis1, dis2, study1, study2,
                            'observed', n_concord, args.qthresh])
    print('Observed concordance... Done.')

    ## Make empirical null for expected number of concordant values
    print('Expected concordance...')
    for n in range(args.nreps):
        print(n),
        # Shuffle the dataframe once per rep
        shuffledeffdir = df.copy().apply(shuffle_col).applymap(sigmap)

        # Then recount the number of concordant directions
        for dis1 in diseases:
            dis1_studies = [i for i in df.columns if i.startswith(dis1)]
            for dis2 in diseases[diseases.index(dis1):]:
                dis2_studies = [i for i in df.columns if i.startswith(dis2)]

                for study1 in dis1_studies:
                    for study2 in dis2_studies:
                        if study1 != study2:
                            n_concord = get_n_concord(
                                shuffledeffdir, study1, study2,
                                args.intersect)
                            dist.append(
                                [dis1, dis2, study1, study2,
                                 'expected', n_concord, args.qthresh])
    print('Expected concordance... Done.')

    # Convert results to dataframe and save output
    distdf = pd.DataFrame(
        data=dist,
        columns=['dis1', 'dis2', 'study1', 'study2',
                 'metric', 'n_concord', 'qthresh'])
    # Make sure edd_singh and noncdi_schubert are correctly labeled
    distdf = distdf\
        .replace('cdi_singh', 'edd_singh')\
        .replace('cdi_schubert2', 'noncdi_schubert')
    distdf.to_csv(args.tidy_fout, sep='\t', index=False)

    pvals = get_pvals_df(distdf)
    pvals.to_csv(args.pvals_fout, sep='\t', index=False)


elif args.method == 'fisher':
    # Fisher exact test, ignoring any genera which are NaN anywhere
    from scipy.stats import fisher_exact

    effdir = np.sign(df).replace(0, np.nan) # anything with zero effect is NaN
    fisherps = []
    for dis1 in diseases:
        dis1_studies = [i for i in df.columns if i.startswith(dis1)]
        for dis2 in diseases[diseases.index(dis1):]:
            dis2_studies = [i for i in df.columns if i.startswith(dis2)]

            for study1 in dis1_studies:
                for study2 in dis2_studies:
                    if study1 != study2:
                        # pd.crosstab automatically drops NaNs
                        mat = pd.crosstab(effdir[study1], effdir[study2])
                        odds, p = fisher_exact(mat)
                        fisherps.append([dis1, dis2, study1, study2, p])
    fisherdf = pd.DataFrame(data=fisherps,
        columns=['dis1', 'dis2', 'study1', 'study2', 'p'])
    fisherdf.to_csv(args.pvals_fout, sep='\t', index=False)

    ## TODO: make redundant?!

elif args.method == 'fisher_all' or args.method == 'cohen':
    from scipy.stats import fisher_exact
    from sklearn.metrics import cohen_kappa_score

    effdir = np.sign(df).replace(0, np.nan) # anything with zero effect is NaN
    ps = []

    for dis1 in diseases:
        dis1_studies = [i for i in df.columns if i.startswith(dis1)]
        for dis2 in diseases[diseases.index(dis1):]:
            dis2_studies = [i for i in df.columns if i.startswith(dis2)]

            for study1 in dis1_studies:
                for study2 in dis2_studies:
                    if study1 != study2:
                        tmp = effdir[[study1, study2]]\
                            .replace(0, np.nan)\
                            .dropna(how='all')
                        #tmp[study1] = \
                        #    np.where(
                        #        tmp[study1].isnull(),
                        #        -1*tmp[study2],
                        #        tmp[study1]
                        #    )
                        #tmp[study2] = \
                        #    np.where(
                        #        tmp[study2].isnull(),
                        #        -1*tmp[study1],
                        #        tmp[study2]
                        #    )
                        if args.method == 'fisher_all':
                            odds, p = fisher_exact(
                                pd.crosstab(tmp[study1], tmp[study2]))
                            ps.append([dis1, dis2, study1, study2, p])
                        elif args.method == 'cohen':
                            k = cohen_kappa_score(tmp[study1], tmp[study2])
                            ps.append([dis1, dis2, study1, study2, k])
    if args.method == 'fisher_all':
        cols = ['dis1', 'dis2', 'study1', 'study2', 'p']
    elif args.method == 'cohen':
        cols = ['dis1', 'dis2', 'study1', 'study2', 'k']

    df = pd.DataFrame(data=ps, columns=cols)
    df.to_csv(args.pvals_fout, sep='\t', index=False)
