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

def get_n_concord(df, col1, col2):
    """
    Get the number of concordant values in col1 and col2 of df.

    Values that are nan in one column and not nan in the other are
    considered discordant. Values that are nan in both columns are
    dropped.

    Note: if you've converted your dataframe of qvalues to +/- 1
    values, this code converts zeros to NaN's so that they are
    not considered concordant (np.nan == np.nan is False,
    but 0 == 0 is True).
    """
    tmp = df[[col1, col2]].dropna(how='all')
    tmp = tmp.replace(0, np.nan)

    # If neither column has anything significant
    if np.all(tmp.isnull()):
        return 0
    else:
        return sum(tmp[col1] == tmp[col2])

p = argparse.ArgumentParser()
p.add_argument('qvals', help='path to file with signed qvalues. Datasets '
    + 'in columns, genera in rows, signed qvalues in values.')
p.add_argument('--qthresh', help='significance threshold for which genera to '
    + 'include in concordance analysis. qthresh = 1 means that effect '
    + 'directions will be compared (regardless of signifiance). qthresh = '
    + '0.05 means that directions of significant genera will be compared. '
    + '[default: %(default)s]', default=1.0, type=float)
p.add_argument('--nreps', help='number of shuffles to build empirical null. '
    + '[default: %(default)s]', default=1000, type=int)
p.add_argument('fout', help='file to write tidy results to.')
args = p.parse_args()

# Read in qvalues
df = pd.read_csv(args.qvals, sep='\t', index_col=0)
# Rename edd_singh to cdi_singh
df = df.rename(columns={'edd_singh': 'cdi_singh',
                        'noncdi_schubert': 'cdi_schubert2'})

# Convert signed qvalues to +/- 1 if significant and 0 if not
sigmap = lambda x: np.sign(x) if abs(x) <= args.qthresh else 0
effdir = df.applymap(sigmap)

# Get diseases to compare between
diseases = list(set([i.split('_')[0] for i in df.columns]))

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
                    n_concord = get_n_concord(effdir, study1, study2)
                    dist.append(
                        [dis1, dis2, study1, study2,
                        'observed', n_concord, args.qthresh])
print('Observed concordance... Done.')

## Make empirical null for expected number of concordant values
print('Expected concordance...')
for n in range(args.nreps):
    print(n),
    # Shuffle the dataframe once per rep
    shuffledeffdir = effdir.copy().apply(shuffle_col)

    # Then recount the number of concordant directions
    for dis1 in diseases:
        dis1_studies = [i for i in df.columns if i.startswith(dis1)]
        for dis2 in diseases[diseases.index(dis1):]:
            dis2_studies = [i for i in df.columns if i.startswith(dis2)]

            for study1 in dis1_studies:
                for study2 in dis2_studies:
                    if study1 != study2:
                        n_concord = get_n_concord(
                            shuffledeffdir, study1, study2)
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
distdf.to_csv(args.fout, sep='\t', index=False)
